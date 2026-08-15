"""GraphistTrainer: trains a GraphistModel and runs downstream latent-space analysis.

Refactored from the original ``VGAEST`` class duplicated across the four
pipeline scripts. Behavior-preserving changes:

- The model is now constructed by the caller (see ``graphist.pipeline``) and
  passed in, rather than being built internally from module-level globals
  (``exp_paths``, ``gmt_paths``, ...) that made ``VGAEST`` impossible to reuse
  across datasets without editing the file.
- ``dec_kl_w`` is dropped: it was accepted as a constructor argument in all
  four original scripts but never actually used anywhere in the loss — an
  unused, misleading parameter.
- The GCN-reconstruction loss weight (``rec_w/100`` in PDAC/BRCA-PACSI vs.
  plain ``rec_w`` in Maynard/BRCA-COMMOT — a copy-paste divergence in the
  originals) is now an explicit ``gcn_rec_w`` argument set per-dataset in
  config, instead of a hardcoded division baked into the loss expression.
- ``Pathway_differential_activity`` (renamed ``pathway_differential_activity``)
  now takes group membership as plain spot-name arrays instead of internally
  re-reading a hardcoded "selected spots" CSV — the caller decides how group
  membership is derived (an annotation column, a Stage-1 Graphist selection
  CSV, etc.).
"""
from typing import Optional, Sequence

import numpy as np
import pandas as pd
import torch
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm

from .data.graph import mask_generator
from .differential import bayesian_differential, differential_activity_report, sample_group_latents
from .losses import gcn_loss, reconstruction_loss
from .model.graphist_model import GraphistModel


class GraphistTrainer:
    """Owns a :class:`GraphistModel` plus the graph tensors it trains against."""

    def __init__(
        self,
        model: GraphistModel,
        X: np.ndarray,
        graph_dict: dict,
        rec_w: float = 1000,
        gcn_w: float = 0.1,
        gcn_rec_w: float = 10,
    ):
        self.model = model
        self.rec_w = rec_w
        self.gcn_w = gcn_w
        self.gcn_rec_w = gcn_rec_w

        self.cell_num = len(X)
        self.X = torch.FloatTensor(np.asarray(X).copy())
        self.input_dim = self.X.shape[1]

        self.adj_norm = graph_dict["adj_norm"]
        self.adj_label = graph_dict["adj_label"]
        self.norm_value = graph_dict["norm_value"]
        self._adj_mask: Optional[torch.Tensor] = None

    def train(self, epochs: int = 200, lr: float = 1e-4, decay: float = 5e-4, n_negatives: int = 1) -> None:
        """Train the model end-to-end (pathway-masked + graph-smoothed + link-reconstruction losses)."""
        optimizer = torch.optim.Adam(self.model.parameters(), lr=lr, weight_decay=decay)
        self.model.train()

        for _ in tqdm(range(epochs)):
            self.model.train()
            optimizer.zero_grad()
            z, mu, logvar, de_feat, x_rec = self.model(self.X, self.adj_norm)

            if self._adj_mask is None:
                self._adj_mask = mask_generator(self.adj_label, self.cell_num, n_negatives=n_negatives)

            loss_gcn = gcn_loss(
                preds=self.model.dc(z, self._adj_mask),
                labels=self._adj_mask.coalesce().values(),
                mu=mu,
                logvar=logvar,
                n_nodes=self.cell_num,
                norm=self.norm_value,
            )
            loss_rec_gcn = reconstruction_loss(de_feat, self.X)
            loss_rec = reconstruction_loss(x_rec, self.X)
            loss = self.rec_w * loss_rec + self.gcn_w * loss_gcn + self.gcn_rec_w * loss_rec_gcn

            loss.backward()
            optimizer.step()

    def process(self):
        """Run a full forward pass in eval mode. Returns (z, de_feat, X_rec) as numpy arrays."""
        self.model.eval()
        z, _, _, de_feat, x_rec = self.model(self.X, self.adj_norm)
        return z.detach().numpy(), de_feat.detach().numpy(), x_rec.detach().numpy()

    def recon(self) -> np.ndarray:
        """Reconstructed expression, standardized (zero mean / unit variance per gene)."""
        self.model.eval()
        _, _, _, _, x_rec = self.model(self.X, self.adj_norm)
        return StandardScaler().fit_transform(x_rec.detach().numpy())

    def pathway_differential_activity(
        self,
        X: pd.DataFrame,
        group1_idx: Sequence[str],
        latent_names: Sequence[str],
        group2_idx: Optional[Sequence[str]] = None,
        group1_name: str = "group1",
        group2_name: str = "rest",
        mode: str = "change",
        delta: float = 2.0,
        fdr_target: float = 0.05,
        n_samples: int = 5000,
        n_permutations: int = 5000,
        use_permutations: bool = True,
        random_seed: Optional[int] = None,
    ) -> pd.DataFrame:
        """Bayesian differential pathway-activity test between two spot groups.

        Parameters
        ----------
        X
            [n_spots, n_genes] expression DataFrame indexed by spot name (same
            index used by ``group1_idx``/``group2_idx``).
        group1_idx
            spot names in group 1.
        group2_idx
            spot names in group 2. If None, group 2 is every spot not in group 1.
        latent_names
            names for each latent dimension (e.g. from
            :func:`graphist.data.pathway_mask.pathway_names`), used as the
            result DataFrame's index.
        """
        with torch.no_grad():
            z = self.model.to_latent(torch.Tensor(X.values), self.adj_norm).numpy()
        z = pd.DataFrame(z, index=X.index.values)

        if group2_idx is None:
            group2_idx = X.index.values[~X.index.isin(group1_idx)]

        z1, z2 = sample_group_latents(z, group1_idx, group2_idx, n_samples=n_samples)
        res = bayesian_differential(
            z1,
            z2,
            mode=mode,
            delta=delta,
            use_permutations=use_permutations,
            n_permutations=n_permutations,
            random_seed=random_seed,
        )
        return differential_activity_report(
            res, latent_names, mode=mode, fdr_target=fdr_target, name_g1=group1_name, name_g2=group2_name
        )
