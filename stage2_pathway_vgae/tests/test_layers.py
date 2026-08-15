import numpy as np
import torch

from graphist.model.layers import CustomizedLinear, GraphConvolution, InnerProductDecoder


def test_customized_linear_forward_shape():
    mask = np.ones((10, 4))
    layer = CustomizedLinear(mask)
    x = torch.randn(5, 10)
    out = layer(x)
    assert out.shape == (5, 4)


def test_customized_linear_masks_forward_weights():
    """Weights outside the mask should be exactly zero even after re-initialization."""
    mask = np.zeros((6, 2))
    mask[:3, 0] = 1  # only first 3 inputs feed output 0
    mask[3:, 1] = 1  # only last 3 inputs feed output 1
    layer = CustomizedLinear(mask)
    effective_weight = layer.weight.data * layer.mask.data
    assert torch.all(effective_weight[0, 3:] == 0)  # output 0 has zero weight from inputs 3-5
    assert torch.all(effective_weight[1, :3] == 0)  # output 1 has zero weight from inputs 0-2


def test_customized_linear_gradient_stays_zero_outside_mask():
    mask = np.zeros((6, 2))
    mask[:3, 0] = 1
    mask[3:, 1] = 1
    layer = CustomizedLinear(mask, bias=False)

    x = torch.randn(4, 6, requires_grad=True)
    out = layer(x)
    out.sum().backward()

    assert torch.all(layer.weight.grad[0, 3:] == 0)
    assert torch.all(layer.weight.grad[1, :3] == 0)
    # weight itself should still be masked to zero after backward (mask isn't updated by grad)
    assert torch.all((layer.weight.data * layer.mask.data)[0, 3:] == 0)


def test_graph_convolution_shape():
    layer = GraphConvolution(in_features=8, out_features=3, dropout=0.0)
    layer.eval()
    x = torch.randn(5, 8)
    adj = torch.eye(5).to_sparse()
    out = layer(x, adj)
    assert out.shape == (5, 3)


def test_inner_product_decoder_self_similarity():
    """sigmoid(z_i . z_i) should be higher than sigmoid(z_i . z_j) for near-orthogonal z_i, z_j."""
    decoder = InnerProductDecoder(dropout=0.0)
    z = torch.tensor([[1.0, 0.0], [0.0, 1.0]])
    # mask selects pairs (0,0) and (0,1)
    indices = torch.tensor([[0, 0], [0, 1]])
    values = torch.tensor([1.0, 1.0])
    mask = torch.sparse_coo_tensor(indices, values, (2, 2))
    result = decoder(z, mask)
    assert result[0] > result[1]
