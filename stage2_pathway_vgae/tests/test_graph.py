import numpy as np
import torch

from graphist.data.graph import generate_adj_mat, graph_construction, mask_generator, preprocess_graph
import scipy.sparse as sp


def test_generate_adj_mat_symmetric_no_self_loops(tiny_adata):
    adj = generate_adj_mat(tiny_adata, include_self=False, n=4)
    assert adj.shape == (40, 40)
    assert np.array_equal(adj, adj.T), "adjacency must be symmetric"
    assert np.all(np.diag(adj) == 0), "no self loops when include_self=False"


def test_generate_adj_mat_include_self(tiny_adata):
    adj = generate_adj_mat(tiny_adata, include_self=True, n=4)
    assert np.all(np.diag(adj) == 1)


def test_generate_adj_mat_each_spot_has_neighbors(tiny_adata):
    adj = generate_adj_mat(tiny_adata, include_self=False, n=4)
    assert np.all(adj.sum(axis=1) >= 4), "every spot should have at least n neighbors after symmetrization"


def test_preprocess_graph_normalization_matches_hand_computation():
    # 3-node path graph: 0-1-2
    adj = sp.coo_matrix(np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]], dtype=float))
    normalized = preprocess_graph(adj).to_dense().numpy()

    adj_plus_i = np.array([[1, 1, 0], [1, 1, 1], [0, 1, 1]], dtype=float)
    deg = adj_plus_i.sum(1)
    d_inv_sqrt = np.diag(deg ** -0.5)
    expected = d_inv_sqrt @ adj_plus_i @ d_inv_sqrt

    np.testing.assert_allclose(normalized, expected, atol=1e-6)


def test_graph_construction_returns_expected_keys_and_shapes(tiny_adata):
    graph_dict = graph_construction(tiny_adata, n=4)
    assert set(graph_dict.keys()) == {"adj_norm", "adj_label", "norm_value"}
    assert graph_dict["adj_norm"].shape == (40, 40)
    assert graph_dict["adj_label"].shape == (40, 40)
    assert graph_dict["norm_value"] > 0


def test_mask_generator_samples_negatives_roughly_matching_positive_count(tiny_adata):
    """mask_generator targets ~1 negative per positive edge (n_negatives=1), but per-node
    sampling (see the docstring quirk in graphist.data.graph.mask_generator) means the
    match isn't exact -- just check it's in the right ballpark, not a strict equality.
    """
    graph_dict = graph_construction(tiny_adata, n=4)
    torch.manual_seed(0)
    mask = mask_generator(graph_dict["adj_label"], cell_num=40, n_negatives=1)
    values = mask.coalesce().values()
    n_pos = int((values == 1).sum())
    n_neg = int((values == 0).sum())
    assert n_pos > 0 and n_neg > 0
    assert 0.5 * n_pos <= n_neg <= 1.5 * n_pos
