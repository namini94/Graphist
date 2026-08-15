"""Bayesian differential pathway-activity testing between two groups of spots.

Adapted from scVI's Bayesian differential expression procedure (Lopez et al.
2018, "vanilla" mode; Boyeau et al. 2019, "change" mode), applied here to
latent pathway-activity scores instead of gene expression. Ported out of
``SVEGA.bayesian_differential``/``STPA.py`` as plain functions operating on
latent arrays directly, so they're testable without a model or graph.
"""
from typing import Tuple

import numpy as np
import pandas as pd


def scale_sampling(arr1: np.ndarray, arr2: np.ndarray, n_perm: int = 1000) -> Tuple[np.ndarray, np.ndarray]:
    """Permutation sampling to better estimate the double integral in the posterior test.

    Inspired by scVI (Lopez et al., 2018): randomly pairs up samples from each
    group (with replacement) to build ``n_perm`` matched comparisons.
    """
    u = np.random.choice(arr1.shape[0], size=n_perm)
    v = np.random.choice(arr2.shape[0], size=n_perm)
    return arr1[u], arr2[v]


def fdr_de_prediction(posterior_probas: np.ndarray, fdr: float = 0.05) -> np.ndarray:
    """Compute posterior expected FDR and flag features as differentially active. From scvi-tools."""
    if posterior_probas.ndim != 1:
        raise ValueError("posterior_probas should be 1-dimensional")
    sorted_genes = np.argsort(-posterior_probas)
    sorted_pgs = posterior_probas[sorted_genes]
    cumulative_fdr = (1.0 - sorted_pgs).cumsum() / (1.0 + np.arange(len(sorted_pgs)))
    d = (cumulative_fdr <= fdr).sum()
    pred_de = np.zeros_like(cumulative_fdr, dtype=bool)
    pred_de[sorted_genes[:d]] = True
    return pred_de


def bayesian_differential(
    z1: np.ndarray,
    z2: np.ndarray,
    mode: str = "change",
    delta: float = 2.0,
    alpha: float = 0.66,
    use_permutations: bool = True,
    n_permutations: int = 5000,
    random_seed: int = None,
) -> dict:
    """Bayesian differential activity test between two groups' latent pathway-activity samples.

    Parameters
    ----------
    z1, z2
        [n_samples, n_latent] latent pathway-activity samples for group 1 / group 2
        (typically drawn with replacement from each group's spots — see
        :func:`sample_group_latents`).
    mode
        ``"vanilla"`` (Lopez et al. 2018) or ``"change"`` (Boyeau et al. 2019).
    delta
        differential-activity threshold used in ``"change"`` mode.
    alpha
        minimum posterior probability (in ``"change"`` mode) to call a pathway
        differentially active alongside the ``delta`` threshold.
    use_permutations
        if True, resample matched pairs via :func:`scale_sampling` first (recommended;
        this is what turns two independent samples into a paired comparison).
    n_permutations
        number of permutation pairs to draw if ``use_permutations``.
    random_seed
        optional seed for reproducibility.

    Returns
    -------
    dict with Bayes factor / mean-difference results, one entry per latent dimension.
    """
    if random_seed is not None:
        np.random.seed(random_seed)
    if mode not in ("vanilla", "change"):
        raise ValueError('mode must be one of "vanilla", "change"')

    epsilon = 1e-12
    if use_permutations:
        z1, z2 = scale_sampling(z1, z2, n_perm=n_permutations)

    if mode == "vanilla":
        p_h1 = np.mean(z1 > z2, axis=0)
        p_h2 = 1.0 - p_h1
        md = np.mean(z1 - z2, axis=0)
        bf = np.log(p_h1 + epsilon) - np.log(p_h2 + epsilon)
        return {"p_h1": p_h1, "p_h2": p_h2, "bayes_factor": bf, "differential_metric": md}

    diffs = z1 - z2
    md = diffs.mean(0)
    p_da = np.mean(np.abs(diffs) > delta, axis=0)
    is_da_alpha = (np.abs(md) > delta) & (p_da > alpha)
    return {
        "p_da": p_da,
        "p_not_da": 1.0 - p_da,
        "bayes_factor": np.log(p_da + epsilon) - np.log((1.0 - p_da) + epsilon),
        f"is_da_alpha_{alpha}": is_da_alpha,
        "differential_metric": md,
        "delta": delta,
    }


def sample_group_latents(
    z: pd.DataFrame, idx1, idx2, n_samples: int = 5000
) -> Tuple[np.ndarray, np.ndarray]:
    """Draw ``n_samples`` spots (with replacement) from each group's latent embedding.

    ``z`` is a [n_spots, n_latent] DataFrame indexed by spot name; ``idx1``/``idx2``
    are the spot names belonging to each group.
    """
    sampled1 = np.random.choice(idx1, n_samples)
    sampled2 = np.random.choice(idx2, n_samples)
    return z.loc[sampled1].values, z.loc[sampled2].values


def differential_activity_report(
    res: dict,
    latent_names,
    mode: str,
    fdr_target: float,
    name_g1: str,
    name_g2: str,
) -> pd.DataFrame:
    """Wrap a :func:`bayesian_differential` result dict into a sorted, annotated DataFrame."""
    df = pd.DataFrame(res, index=latent_names)
    sort_key = "p_da" if mode == "change" else "bayes_factor"
    df = df.sort_values(by=sort_key, ascending=False)
    if mode == "change":
        df[f"is_da_fdr_{fdr_target}"] = fdr_de_prediction(df["p_da"].values, fdr=fdr_target)
    df["comparison"] = f"{name_g1} vs. {name_g2}"
    df["group1"] = name_g1
    df["group2"] = name_g2
    return df
