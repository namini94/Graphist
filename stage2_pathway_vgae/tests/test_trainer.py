import numpy as np
import torch

from graphist.data.graph import graph_construction
from graphist.data.pathway_mask import create_mask
from graphist.model.graphist_model import GraphistModel
from graphist.trainer import GraphistTrainer
from graphist.utils import seed_everything


def _build_trainer(tiny_adata, tiny_gmt_path):
    seed_everything(0)
    graph_dict = graph_construction(tiny_adata, n=4)
    features = tiny_adata.var_names.tolist()
    mask = create_mask(features, tiny_gmt_path, add_nodes=1)
    model = GraphistModel(input_dim=len(features), n_gmvs=mask.shape[1], mask=mask, gcn_hidden1=16)
    trainer = GraphistTrainer(model, np.asarray(tiny_adata.X), graph_dict, rec_w=10, gcn_w=0.1, gcn_rec_w=1)
    return trainer


def test_train_runs_without_nan(tiny_adata, tiny_gmt_path):
    trainer = _build_trainer(tiny_adata, tiny_gmt_path)
    trainer.train(epochs=2, lr=1e-3)
    z, de_feat, x_rec = trainer.process()
    assert not np.isnan(z).any()
    assert not np.isnan(de_feat).any()
    assert not np.isnan(x_rec).any()


def test_process_output_shapes(tiny_adata, tiny_gmt_path):
    trainer = _build_trainer(tiny_adata, tiny_gmt_path)
    trainer.train(epochs=1, lr=1e-3)
    z, de_feat, x_rec = trainer.process()
    n_spots, n_genes = tiny_adata.shape
    assert z.shape[0] == n_spots
    assert de_feat.shape == (n_spots, n_genes)
    assert x_rec.shape == (n_spots, n_genes)


def test_recon_is_standardized(tiny_adata, tiny_gmt_path):
    trainer = _build_trainer(tiny_adata, tiny_gmt_path)
    trainer.train(epochs=1, lr=1e-3)
    out = trainer.recon()
    # StandardScaler -> each gene column should have ~zero mean
    assert np.allclose(out.mean(axis=0), 0, atol=1e-6)


def test_training_reduces_reconstruction_loss(tiny_adata, tiny_gmt_path):
    """Sanity check that the loss trends down over a short training run."""
    trainer = _build_trainer(tiny_adata, tiny_gmt_path)
    trainer.train(epochs=1, lr=1e-2)
    _, _, x_rec_before = trainer.process()
    loss_before = np.mean((x_rec_before - np.asarray(tiny_adata.X)) ** 2)

    trainer.train(epochs=30, lr=1e-2)
    _, _, x_rec_after = trainer.process()
    loss_after = np.mean((x_rec_after - np.asarray(tiny_adata.X)) ** 2)

    assert loss_after < loss_before
