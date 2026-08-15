"""GraphistModel: the pathway-activity variational graph autoencoder.

This is the refactored, renamed ``SVEGA`` class from the four original
per-dataset scripts. Behavior-preserving changes:

- Takes an already-computed ``mask`` array instead of reading an expression
  CSV off disk itself (the original constructor called
  ``pd.read_csv(exp_paths)`` internally purely to recover the gene name list
  for :func:`~graphist.data.pathway_mask.create_mask` — that's now the
  caller's responsibility, via :mod:`graphist.pipeline`).
- Drops ``differential_activity``/the single-``bayesian_differential``-without-``adj``
  code path: that method was hardcoded against an unrelated single-cell
  drug-perturbation dataset (``metadata_Belin.csv``) inherited from the
  original single-cell VEGA project and was never called by any of the four
  spatial pipeline scripts. The real, used differential-activity logic lives
  in :mod:`graphist.differential` (pure functions) and is orchestrated by
  :class:`~graphist.trainer.GraphistTrainer`.
"""
import os

import numpy as np
import torch
import torch.nn.functional as F
from torch import nn

from .decoder import PathwayMaskedDecoder
from .layers import GraphConvolution, InnerProductDecoder


class GraphistModel(nn.Module):
    """GCN encoder + three decoders: pathway-masked, graph-smoothed, and link-reconstruction.

    Parameters
    ----------
    input_dim
        number of (HVG) genes.
    n_gmvs
        number of latent pathway-activity dimensions (n_pathways + add_nodes).
    mask
        [n_genes, n_gmvs] gene x pathway membership mask (see
        :func:`graphist.data.pathway_mask.create_mask`).
    dropout / gcn_hidden1 / p_drop
        encoder hyperparameters.
    positive_decoder
        constrain the pathway decoder's weights to be non-negative.
    """

    def __init__(
        self,
        input_dim: int,
        n_gmvs: int,
        mask: np.ndarray,
        dropout: float = 0.1,
        positive_decoder: bool = True,
        gcn_hidden1: int = 800,
        p_drop: float = 0.2,
    ):
        super().__init__()

        # NOTE: construction order matches the original SVEGA class exactly
        # (gc_decoder, then gc1/gc2/gc3/dc, then decoder last). Each of these
        # submodules draws from the global torch RNG during weight init
        # (xavier_uniform_ / uniform_), so under a fixed seed, construction
        # order determines which random draws each layer gets -- reordering
        # would silently change the trained model despite an identical seed.
        self.gc_decoder = GraphConvolution(n_gmvs, input_dim, p_drop, act=lambda x: x)

        self.gc1 = GraphConvolution(input_dim, gcn_hidden1, p_drop, act=F.relu)
        self.gc2 = GraphConvolution(gcn_hidden1, n_gmvs, p_drop, act=lambda x: x)
        self.gc3 = GraphConvolution(gcn_hidden1, n_gmvs, p_drop, act=lambda x: x)
        self.dc = InnerProductDecoder(p_drop, act=lambda x: x)

        # Pathway-masked expression reconstruction decoder
        self.decoder = PathwayMaskedDecoder(mask=mask.T, positive_decoder=positive_decoder)

    def encode(self, X: torch.Tensor, adj: torch.Tensor):
        hidden1 = self.gc1(X, adj)
        return self.gc2(hidden1, adj), self.gc3(hidden1, adj)

    def sample_latent(self, mu: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
        std = logvar.mul(0.5).exp_()
        eps = torch.FloatTensor(std.size()).normal_()
        return eps.mul_(std).add_(mu)

    @torch.no_grad()
    def to_latent(self, X: torch.Tensor, adj: torch.Tensor) -> torch.Tensor:
        """Same as :meth:`encode` + :meth:`sample_latent`, but returns only ``z``."""
        mu, logvar = self.encode(X, adj)
        return self.sample_latent(mu, logvar)

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        return self.decoder(z)

    def forward(self, x: torch.Tensor, adj: torch.Tensor):
        mu, logvar = self.encode(x, adj)
        z = self.sample_latent(mu, logvar)
        de_feat = self.gc_decoder(z, adj)
        x_rec = self.decode(z)
        return z, mu, logvar, de_feat, x_rec

    def save(self, path: str, overwrite: bool = True) -> None:
        """Save model parameters to ``path`` (creates the directory if needed)."""
        if not os.path.exists(path) or overwrite:
            os.makedirs(path, exist_ok=overwrite)
        else:
            raise ValueError(f"{path} already exists. Please provide a non-existing directory for saving.")
        torch.save(self.state_dict(), os.path.join(path, "graphist_params.pt"))
        print(f"Model files saved at {path}")

    def load(self, path: str) -> None:
        """Load model parameters previously written by :meth:`save`."""
        state_dict = torch.load(os.path.join(path, "graphist_params.pt"))
        self.load_state_dict(state_dict)
        print(f"Model loaded from {path}")
