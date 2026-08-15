"""Spatial neighbor graph construction for the GCN encoder/decoder.

Ported from the shared boilerplate duplicated across the four original
per-dataset pipeline scripts (generate_adj_mat / preprocess_graph /
graph_construction / mask_generator were byte-identical or near-identical in
all four). Two genuinely dead helpers from the originals were dropped:
``graph_computing`` (an alternate KNN implementation, never called — every
script always used ``generate_adj_mat`` via ``mode='KNN'``) and
``combine_graph_dict`` (multi-sample graph merging, never called by any of
the single-sample per-dataset scripts).
"""
from typing import Optional

import numpy as np
import scipy.sparse as sp
import torch
from anndata import AnnData
from sklearn import metrics


def generate_adj_mat(adata: AnnData, include_self: bool = False, n: int = 6) -> np.ndarray:
    """Build a symmetric binary KNN adjacency matrix from ``adata.obsm['spatial']``."""
    assert "spatial" in adata.obsm, "AnnData object should provide spatial information"

    dist = metrics.pairwise_distances(adata.obsm["spatial"])

    adj_mat = np.zeros((len(adata), len(adata)))
    for i in range(len(adata)):
        n_neighbors = np.argsort(dist[i, :])[: n + 1]
        adj_mat[i, n_neighbors] = 1

    if not include_self:
        x, y = np.diag_indices_from(adj_mat)
        adj_mat[x, y] = 0

    adj_mat = adj_mat + adj_mat.T
    adj_mat = adj_mat > 0
    adj_mat = adj_mat.astype(np.int64)

    return adj_mat


def generate_adj_mat_radius(adata: AnnData, max_dist: float) -> np.ndarray:
    """Build a binary adjacency matrix by thresholding pairwise Euclidean distance."""
    assert "spatial" in adata.obsm, "AnnData object should provide spatial information"

    dist = metrics.pairwise_distances(adata.obsm["spatial"], metric="euclidean")
    adj_mat = dist < max_dist
    adj_mat = adj_mat.astype(np.int64)
    return adj_mat


def sparse_mx_to_torch_sparse_tensor(sparse_mx: sp.spmatrix) -> torch.Tensor:
    """Convert a scipy sparse matrix to a torch sparse tensor."""
    sparse_mx = sparse_mx.tocoo().astype(np.float32)
    indices = torch.from_numpy(np.vstack((sparse_mx.row, sparse_mx.col)).astype(np.int64))
    values = torch.from_numpy(sparse_mx.data)
    shape = torch.Size(sparse_mx.shape)
    return torch.sparse_coo_tensor(indices, values, shape)


def preprocess_graph(adj: sp.spmatrix) -> torch.Tensor:
    """Symmetric degree normalization: D^-1/2 (A + I) D^-1/2."""
    adj_ = adj + sp.eye(adj.shape[0])
    rowsum = np.array(adj_.sum(1))
    degree_mat_inv_sqrt = sp.diags(np.power(rowsum, -0.5).flatten())
    adj_normalized = adj_.dot(degree_mat_inv_sqrt).transpose().dot(degree_mat_inv_sqrt).tocoo()
    return sparse_mx_to_torch_sparse_tensor(adj_normalized)


def mask_generator(adj_label: torch.Tensor, cell_num: int, n_negatives: int = 1) -> torch.Tensor:
    """Build a (positive-edge + random-negative-edge) mask for the VGAE link-reconstruction loss.

    ``cell_num`` is passed explicitly (rather than inferred) since this is also called
    once per training step on the full graph inside ``GraphistTrainer.train``.

    Known quirk, preserved intentionally
    -------------------------------------
    For each node, this samples ``len(neighbors) * n_negatives`` "negative" target
    node ids from ``torch.randperm(len(non_neighbor))`` — but then uses those
    permutation *positions* directly as target node ids, rather than indexing
    back into ``non_neighbor`` to recover the actual non-neighbor node ids
    (i.e. it never does ``non_neighbor[indices[:n_selected]]``). In a sparse
    graph this still lands close to a uniform sample of node ids, but it is not
    a strict guarantee that every "negative" edge is a true non-edge, and the
    resulting negative-edge count is not guaranteed to exactly equal the
    positive-edge count (some nodes can run out of distinct sampled positions).
    This exists identically in all four original per-dataset scripts (see
    ``legacy/*.py``) and is preserved here bug-for-bug so that a refactored
    training run reproduces the original models' loss dynamics; it is *not*
    silently "fixed", since that would change what gets learned.
    """
    idx = adj_label.indices()

    list_non_neighbor = []
    for i in range(cell_num):
        neighbor = idx[1, torch.where(idx[0, :] == i)[0]]
        n_selected = len(neighbor) * n_negatives

        total_idx = torch.arange(0, cell_num, dtype=torch.float32)
        non_neighbor = total_idx[~torch.isin(total_idx, neighbor)]
        indices = torch.randperm(len(non_neighbor), dtype=torch.float32)
        random_non_neighbor = indices[:n_selected]
        list_non_neighbor.append(random_non_neighbor)

    x = torch.repeat_interleave(idx[0], n_negatives)
    y = torch.concat(list_non_neighbor)

    neg_indices = torch.stack([x, y])
    all_indices = torch.concat([idx, neg_indices], axis=1)

    values = torch.concat([adj_label.values(), torch.zeros(len(x), dtype=torch.float32)])
    adj_mask = torch.sparse_coo_tensor(all_indices, values)

    return adj_mask


def graph_construction(adata: AnnData, n: int = 6, dmax: float = 50, mode: str = "KNN") -> dict:
    """Build the normalized/label adjacency tensors consumed by GraphistTrainer.

    Returns a dict with:
        adj_norm: symmetrically-normalized adjacency (for the GCN encoder/decoders)
        adj_label: adjacency + self loops, as the link-reconstruction target
        norm_value: class-imbalance correction weight for the reconstruction BCE loss
    """
    if mode == "KNN":
        adj_m1 = generate_adj_mat(adata, include_self=False, n=n)
    else:
        adj_m1 = generate_adj_mat_radius(adata, dmax)
    adj_m1 = sp.coo_matrix(adj_m1)

    # Store original adjacency matrix (without diagonal entries) for later
    adj_m1 = adj_m1 - sp.dia_matrix((adj_m1.diagonal()[np.newaxis, :], [0]), shape=adj_m1.shape)
    adj_m1.eliminate_zeros()

    adj_norm_m1 = preprocess_graph(adj_m1)
    adj_m1 = adj_m1 + sp.eye(adj_m1.shape[0])

    adj_m1 = adj_m1.tocoo()
    shape = adj_m1.shape
    values = adj_m1.data
    indices = np.stack([adj_m1.row, adj_m1.col])
    adj_label_m1 = torch.sparse_coo_tensor(indices, values, shape)

    norm_m1 = adj_m1.shape[0] * adj_m1.shape[0] / float(
        (adj_m1.shape[0] * adj_m1.shape[0] - adj_m1.sum()) * 2
    )

    return {
        "adj_norm": adj_norm_m1,
        "adj_label": adj_label_m1.coalesce(),
        "norm_value": norm_m1,
    }
