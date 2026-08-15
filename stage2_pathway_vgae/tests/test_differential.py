import numpy as np
import pandas as pd

from graphist.differential import (
    bayesian_differential,
    differential_activity_report,
    fdr_de_prediction,
    sample_group_latents,
    scale_sampling,
)


def test_bayesian_differential_vanilla_recovers_shifted_direction():
    rng = np.random.default_rng(0)
    z1 = rng.normal(loc=5.0, scale=0.5, size=(200, 3))  # dim 0 shifted up in group 1
    z2 = rng.normal(loc=0.0, scale=0.5, size=(200, 3))
    res = bayesian_differential(z1, z2, mode="vanilla", use_permutations=False, random_seed=0)
    assert res["bayes_factor"][0] > 0
    assert res["differential_metric"][0] > 4  # roughly the 5.0 - 0.0 shift


def test_bayesian_differential_change_mode_flags_shifted_dim():
    rng = np.random.default_rng(0)
    # Only dimension 0 differs between groups; dims 1-2 are drawn from the same
    # distribution in both groups and should NOT be flagged as differentially active.
    z1 = rng.normal(loc=0.0, scale=0.5, size=(200, 3))
    z1[:, 0] += 5.0
    z2 = rng.normal(loc=0.0, scale=0.5, size=(200, 3))
    res = bayesian_differential(z1, z2, mode="change", delta=1.0, use_permutations=False, random_seed=0)
    assert res["is_da_alpha_0.66"][0]  # shifted dimension flagged differentially active
    assert not res["is_da_alpha_0.66"][1]  # unshifted dims should not be flagged
    assert not res["is_da_alpha_0.66"][2]


def test_bayesian_differential_symmetric_groups_near_zero_effect():
    rng = np.random.default_rng(1)
    z1 = rng.normal(loc=0.0, scale=1.0, size=(500, 2))
    z2 = rng.normal(loc=0.0, scale=1.0, size=(500, 2))
    res = bayesian_differential(z1, z2, mode="vanilla", use_permutations=False, random_seed=1)
    assert np.abs(res["differential_metric"]).max() < 0.5


def test_scale_sampling_output_shapes():
    arr1 = np.arange(20).reshape(10, 2)
    arr2 = np.arange(20, 40).reshape(10, 2)
    s1, s2 = scale_sampling(arr1, arr2, n_perm=100)
    assert s1.shape == (100, 2)
    assert s2.shape == (100, 2)


def test_fdr_de_prediction_thresholding():
    # Strong evidence for the first 3, weak/none for the rest
    probs = np.array([0.99, 0.98, 0.97, 0.4, 0.3, 0.1])
    flagged = fdr_de_prediction(probs, fdr=0.05)
    assert flagged[:3].all()
    assert not flagged[3:].any()


def test_sample_group_latents_uses_correct_groups():
    z = pd.DataFrame({"lat0": [1, 2, 3, 100, 101, 102]}, index=["a", "b", "c", "d", "e", "f"])
    z1, z2 = sample_group_latents(z, idx1=["a", "b", "c"], idx2=["d", "e", "f"], n_samples=50)
    assert z1.max() < 10
    assert z2.min() > 50


def test_differential_activity_report_is_sorted_and_annotated():
    res = {"bayes_factor": np.array([0.1, 5.0, 2.0]), "differential_metric": np.array([0.1, 3.0, 1.0])}
    df = differential_activity_report(
        res, latent_names=["P1", "P2", "P3"], mode="vanilla", fdr_target=0.05, name_g1="g1", name_g2="g2"
    )
    assert df.index[0] == "P2"  # highest bayes_factor sorted first
    assert (df["comparison"] == "g1 vs. g2").all()
