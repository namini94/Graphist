"""Pathway-masked decoder: reconstructs gene expression from pathway-activity latents.

Adapted from VEGA's ``LinearDecoderSCVI``-inspired masked decoder
(https://github.com/YosefLab/scvi-tools) and the "Customized Linear" masked
layer (https://github.com/uchida-takumi/CustomizedLinear). The original had a
module-level ``n_out = 5000`` global feeding ``SparseLayer``'s optional
BatchNorm dimension, hardcoded independently of the mask actually passed in;
that's now derived from the mask shape (see ``SparseLayer`` in ``layers.py``),
so this module no longer needs to know the gene count in advance.
"""
from typing import Iterable, Optional

import numpy as np
import torch
from torch import nn

from .layers import CustomizedLinear, SparseLayer


class PathwayMaskedDecoder(nn.Module):
    """Decoder for pathway-activity latents (log/z-scored expression data).

    Parameters
    ----------
    mask
        [n_latent, n_genes] gene-module membership matrix (transposed relative
        to the [n_genes, n_latent] mask used for the encoder side).
    n_cat_list
        list encoding number of categories for each covariate.
    positive_decoder
        whether to constrain decoder weights to non-negative values, so that
        pathway "activity" has a consistent sign (higher = more active).
    """

    def __init__(
        self,
        mask: np.ndarray,
        n_cat_list: Optional[Iterable[int]] = None,
        positive_decoder: bool = True,
    ):
        super().__init__()
        self.n_input = mask.shape[0]
        self.n_output = mask.shape[1]
        self.decoder = SparseLayer(
            mask,
            n_cat_list=n_cat_list,
            use_batch_norm=False,
            use_layer_norm=False,
            bias=True,
            dropout_rate=0,
        )
        if positive_decoder:
            self._positive_weights()

    def forward(self, x: torch.Tensor, *cat_list: int) -> torch.Tensor:
        return self.decoder(x, *cat_list)

    def _get_weights(self) -> torch.Tensor:
        """Weight matrix of the masked linear decoder (for inspection/regularization)."""
        return self.decoder.sparse_layer[0].weight

    def _positive_weights(self, use_softplus: bool = False) -> None:
        """Clamp (or softplus-transform) decoder weights to be non-negative."""
        w = self._get_weights()
        if use_softplus:
            w.data = nn.functional.softplus(w.data)
        else:
            w.data = w.data.clamp(0)
