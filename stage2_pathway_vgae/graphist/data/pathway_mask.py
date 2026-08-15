"""Gene x pathway mask construction from .gmt pathway files.

Ported from STPA.py / the duplicated block in each of the four original
pipeline scripts. The only behavior change: ``create_mask`` no longer writes
a debug CSV to a hardcoded ``/Users/naminiyakan/...`` path as a side effect.
"""
from collections import OrderedDict
from typing import List, Union

import numpy as np


def read_gmt(fname: str, sep: str = "\t", min_g: int = 0, max_g: int = 5000) -> "OrderedDict[str, list]":
    """Read a .gmt file into an ordered dict of pathway_name -> [gene, gene, ...].

    min_g/max_g optionally filter out pathways outside a gene-set size range.
    """
    dict_gmv: "OrderedDict[str, list]" = OrderedDict()
    with open(fname) as f:
        for line in f.readlines():
            line = line.strip()
            val = line.split(sep)
            if min_g <= len(val[2:]) <= max_g:
                dict_gmv[val[0]] = val[2:]
    return dict_gmv


def write_gmt(dict_obj: dict, path_gmt: str, sep: str = "\t", second_col: bool = True) -> None:
    """Write a pathway_name -> [genes] dict out to .gmt format."""
    with open(path_gmt, "w") as f:
        for k, v in dict_obj.items():
            if second_col:
                to_write = sep.join([k, "SECOND_COL"] + v) + "\n"
            else:
                to_write = sep.join([k] + v) + "\n"
            f.write(to_write)


def make_pathway_mask(feature_list: List[str], dict_gmv: "OrderedDict[str, list]", add_nodes: int) -> np.ndarray:
    """Build a [n_genes, n_pathways + add_nodes] binary membership mask.

    (i, j) = 1 if gene i belongs to pathway j, 0 otherwise. The last
    ``add_nodes`` columns are fully-connected (all ones) "unannotated" latent
    nodes that can absorb residual variance not captured by any pathway.
    """
    assert isinstance(dict_gmv, OrderedDict), "dict_gmv must be an OrderedDict so column order is stable"
    p_mask = np.zeros((len(feature_list), len(dict_gmv)))
    for j, k in enumerate(dict_gmv.keys()):
        genes_in_pathway = set(dict_gmv[k])
        for i, gene in enumerate(feature_list):
            if gene in genes_in_pathway:
                p_mask[i, j] = 1.0

    unannotated = np.ones((p_mask.shape[0], add_nodes))
    p_mask = np.hstack((p_mask, unannotated))
    return p_mask


def create_mask(
    feature_list: List[str],
    gmt_paths: Union[str, List[str]],
    add_nodes: int = 1,
    min_genes: int = 0,
    max_genes: int = 1000,
) -> np.ndarray:
    """Build the gene x pathway mask used by :class:`~graphist.model.decoder.PathwayMaskedDecoder`.

    Parameters
    ----------
    feature_list
        Gene names in the (HVG-filtered) dataset, in column order.
    gmt_paths
        One or several .gmt pathway files.
    add_nodes
        Number of additional, fully-connected "unannotated" latent nodes.
    min_genes / max_genes
        Gene-set size filters applied per pathway.
    """
    if isinstance(gmt_paths, str):
        gmt_paths = [gmt_paths]

    dict_gmv: "OrderedDict[str, list]" = OrderedDict()
    for f in gmt_paths:
        dict_gmv.update(read_gmt(f, sep="\t", min_g=min_genes, max_g=max_genes))

    return make_pathway_mask(feature_list=feature_list, dict_gmv=dict_gmv, add_nodes=add_nodes)


def pathway_names(gmt_paths: Union[str, List[str]], add_nodes: int = 1, min_genes: int = 0, max_genes: int = 1000) -> List[str]:
    """Return latent-dimension names: pathway names followed by UNANNOTATED_0..add_nodes-1.

    Column order must match :func:`create_mask` exactly since both are derived
    from the same ``read_gmt`` call over the same gmt file(s).
    """
    if isinstance(gmt_paths, str):
        gmt_paths = [gmt_paths]

    dict_gmv: "OrderedDict[str, list]" = OrderedDict()
    for f in gmt_paths:
        dict_gmv.update(read_gmt(f, sep="\t", min_g=min_genes, max_g=max_genes))

    return list(dict_gmv.keys()) + [f"UNANNOTATED_{k}" for k in range(add_nodes)]
