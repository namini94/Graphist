import numpy as np
import torch

from graphist.data.graph import graph_construction
from graphist.data.pathway_mask import create_mask
from graphist.model.graphist_model import GraphistModel


def _build_model_and_inputs(tiny_adata, tiny_gmt_path, add_nodes=1):
    graph_dict = graph_construction(tiny_adata, n=4)
    features = tiny_adata.var_names.tolist()
    mask = create_mask(features, tiny_gmt_path, add_nodes=add_nodes)
    n_gmvs = mask.shape[1]
    model = GraphistModel(input_dim=len(features), n_gmvs=n_gmvs, mask=mask, gcn_hidden1=16)
    X = torch.FloatTensor(np.asarray(tiny_adata.X))
    return model, X, graph_dict["adj_norm"], n_gmvs


def test_forward_pass_shapes(tiny_adata, tiny_gmt_path):
    model, X, adj_norm, n_gmvs = _build_model_and_inputs(tiny_adata, tiny_gmt_path)
    z, mu, logvar, de_feat, x_rec = model(X, adj_norm)
    n_spots, n_genes = X.shape
    assert z.shape == (n_spots, n_gmvs)
    assert mu.shape == (n_spots, n_gmvs)
    assert logvar.shape == (n_spots, n_gmvs)
    assert de_feat.shape == (n_spots, n_genes)
    assert x_rec.shape == (n_spots, n_genes)


def test_to_latent_matches_forward_dimension(tiny_adata, tiny_gmt_path):
    model, X, adj_norm, n_gmvs = _build_model_and_inputs(tiny_adata, tiny_gmt_path)
    z = model.to_latent(X, adj_norm)
    assert z.shape == (X.shape[0], n_gmvs)


def test_positive_decoder_clamps_weights_nonnegative(tiny_adata, tiny_gmt_path):
    model, _, _, _ = _build_model_and_inputs(tiny_adata, tiny_gmt_path)
    w = model.decoder._get_weights()
    assert torch.all(w.data >= 0)


def test_save_load_round_trip(tiny_adata, tiny_gmt_path, tmp_path):
    model, X, adj_norm, n_gmvs = _build_model_and_inputs(tiny_adata, tiny_gmt_path)
    save_dir = str(tmp_path / "model_out")
    model.save(save_dir)

    reloaded = GraphistModel(
        input_dim=X.shape[1],
        n_gmvs=n_gmvs,
        mask=create_mask(tiny_adata.var_names.tolist(), tiny_gmt_path, add_nodes=1),
        gcn_hidden1=16,
    )
    reloaded.load(save_dir)

    for p1, p2 in zip(model.parameters(), reloaded.parameters()):
        assert torch.allclose(p1, p2)
