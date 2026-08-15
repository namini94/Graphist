"""Core neural network layers: GCN convolution, inner-product graph decoder,
and the mask-constrained ("customized") linear layer used by the pathway decoder.
"""
import math

import numpy as np
import torch
import torch.nn.functional as F
from torch import nn
from torch.nn.modules.module import Module
from torch.nn.parameter import Parameter


def one_hot(index: torch.Tensor, n_cat: int) -> torch.Tensor:
    """Minimal one-hot encoding helper for categorical covariates.

    Replaces a dependency on ``scvi.nn.one_hot`` (dropped so graphist no longer
    requires installing scvi-tools). Only exercised if a caller actually passes
    categorical covariates to :class:`SparseLayer` — none of the built-in
    dataset configs do.
    """
    onehot = torch.zeros(index.size(0), n_cat, device=index.device)
    onehot.scatter_(1, index.long(), 1)
    return onehot.type(torch.float32)


class GraphConvolution(Module):
    """Simple GCN layer, similar to https://arxiv.org/abs/1609.02907."""

    def __init__(self, in_features: int, out_features: int, dropout: float = 0.0, act=F.relu):
        super().__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.dropout = dropout
        self.act = act
        self.weight = Parameter(torch.FloatTensor(in_features, out_features))
        self.reset_parameters()

    def reset_parameters(self):
        nn.init.xavier_uniform_(self.weight)

    def forward(self, input: torch.Tensor, adj: torch.Tensor) -> torch.Tensor:
        input = F.dropout(input, self.dropout, self.training)
        support = torch.mm(input, self.weight)
        output = torch.spmm(adj, support)
        return self.act(output)


class InnerProductDecoder(nn.Module):
    """Reconstructs graph edges from latent embeddings via a (masked) inner product."""

    def __init__(self, dropout: float, act=torch.sigmoid):
        super().__init__()
        self.dropout = dropout
        self.act = act

    def forward(self, z: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
        col = mask.coalesce().indices()[0]
        row = mask.coalesce().indices()[1]
        return self.act(torch.sum(z[col] * z[row], axis=1))


class CustomizedLinearFunction(torch.autograd.Function):
    """Autograd function for a linear layer whose weights are masked to a fixed sparsity pattern."""

    @staticmethod
    def forward(ctx, input, weight, bias=None, mask=None):
        if mask is not None:
            weight = weight * mask
        output = input.mm(weight.t())
        if bias is not None:
            output += bias.unsqueeze(0).expand_as(output)
        ctx.save_for_backward(input, weight, bias, mask)
        return output

    @staticmethod
    def backward(ctx, grad_output):
        input, weight, bias, mask = ctx.saved_tensors
        grad_input = grad_weight = grad_bias = grad_mask = None

        if ctx.needs_input_grad[0]:
            grad_input = grad_output.mm(weight)
        if ctx.needs_input_grad[1]:
            grad_weight = grad_output.t().mm(input)
            if mask is not None:
                grad_weight = grad_weight * mask
        if ctx.needs_input_grad[2]:
            grad_bias = grad_output.sum(0).squeeze(0)

        return grad_input, grad_weight, grad_bias, grad_mask


class CustomizedLinear(nn.Module):
    """Linear layer constrained to a fixed gene x pathway connectivity mask.

    Adapted from https://github.com/uchida-takumi/CustomizedLinear (as in the
    original VEGA-derived code) — weights outside ``mask`` are always zero,
    both in the forward pass and in the accumulated gradient.
    """

    def __init__(self, mask, bias: bool = True):
        super().__init__()
        self.input_features = mask.shape[0]
        self.output_features = mask.shape[1]
        if isinstance(mask, torch.Tensor):
            self.mask = mask.type(torch.float).t()
        else:
            self.mask = torch.tensor(mask, dtype=torch.float).t()
        self.mask = nn.Parameter(self.mask, requires_grad=False)

        self.weight = nn.Parameter(torch.Tensor(self.output_features, self.input_features))
        if bias:
            self.bias = nn.Parameter(torch.Tensor(self.output_features))
        else:
            self.register_parameter("bias", None)
        self.reset_parameters()

        self.weight.data = self.weight.data * self.mask

    def reset_parameters(self):
        stdv = 1.0 / math.sqrt(self.weight.size(1))
        self.weight.data.uniform_(-stdv, stdv)
        if self.bias is not None:
            self.bias.data.uniform_(-stdv, stdv)

    def reset_params_pos(self):
        """Same as reset_parameters, but only initializes to positive values."""
        stdv = 1.0 / math.sqrt(self.weight.size(1))
        self.weight.data.uniform_(0, stdv)
        if self.bias is not None:
            self.bias.data.uniform_(-stdv, stdv)

    def forward(self, input: torch.Tensor) -> torch.Tensor:
        return CustomizedLinearFunction.apply(input, self.weight, self.bias, self.mask)

    def extra_repr(self) -> str:
        return "input_features={}, output_features={}, bias={}".format(
            self.input_features, self.output_features, self.bias is not None
        )


class SparseLayer(nn.Module):
    """Single masked linear layer with optional norm/activation/dropout and categorical covariates.

    Inspired by scvi-tools' ``FCLayers`` but constrained to exactly one layer.
    """

    def __init__(
        self,
        mask,
        n_cat_list=None,
        n_continuous_cov: int = 0,
        use_activation: bool = False,
        use_batch_norm: bool = False,
        use_layer_norm: bool = False,
        bias: bool = True,
        dropout_rate: float = 0.1,
        activation_fn=None,
    ):
        super().__init__()
        if n_cat_list is not None:
            self.n_cat_list = [n_cat if n_cat > 1 else 0 for n_cat in n_cat_list]
        else:
            self.n_cat_list = []
        self.n_continuous_cov = n_continuous_cov
        self.cat_dim = sum(self.n_cat_list)
        n_out = mask.shape[1]
        mask_with_cov = np.vstack((mask, np.ones((self.n_continuous_cov + self.cat_dim, n_out))))
        self.sparse_layer = nn.Sequential(
            CustomizedLinear(mask_with_cov),
            nn.BatchNorm1d(n_out, momentum=0.01, eps=0.001) if use_batch_norm else None,
            nn.LayerNorm(n_out, elementwise_affine=False) if use_layer_norm else None,
            activation_fn() if use_activation else None,
            nn.Dropout(p=dropout_rate) if dropout_rate > 0 else None,
        )

    def forward(self, x: torch.Tensor, *cat_list: int) -> torch.Tensor:
        one_hot_cat_list = []
        if len(self.n_cat_list) > len(cat_list):
            raise ValueError("nb. categorical args provided doesn't match init. params.")
        for n_cat, cat in zip(self.n_cat_list, cat_list):
            if n_cat and cat is None:
                raise ValueError("cat not provided while n_cat != 0 in init. params.")
            if n_cat > 1:
                one_hot_cat = one_hot(cat, n_cat) if cat.size(1) != n_cat else cat
                one_hot_cat_list += [one_hot_cat]

        for layer in self.sparse_layer:
            if layer is None:
                continue
            if isinstance(layer, nn.BatchNorm1d):
                if x.dim() == 3:
                    x = torch.cat([(layer(slice_x)).unsqueeze(0) for slice_x in x], dim=0)
                else:
                    x = layer(x)
            else:
                if isinstance(layer, CustomizedLinear):
                    if x.dim() == 3:
                        one_hot_cat_list_layer = [
                            o.unsqueeze(0).expand((x.size(0), o.size(0), o.size(1)))
                            for o in one_hot_cat_list
                        ]
                    else:
                        one_hot_cat_list_layer = one_hot_cat_list
                    x = torch.cat((x, *one_hot_cat_list_layer), dim=-1)
                x = layer(x)
        return x
