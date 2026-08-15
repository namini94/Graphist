import torch

from graphist.losses import gcn_loss, reconstruction_loss


def test_reconstruction_loss_zero_when_identical():
    x = torch.randn(5, 4)
    assert reconstruction_loss(x, x).item() == 0.0


def test_reconstruction_loss_positive_when_different():
    x = torch.randn(5, 4)
    y = x + 1.0
    assert reconstruction_loss(y, x).item() > 0.0


def test_gcn_loss_kld_zero_at_standard_normal():
    """mu=0, logvar=0 (i.e. std=1) exactly matches the N(0,1) prior -> KL term is 0."""
    n_nodes = 5
    mu = torch.zeros(n_nodes, 3)
    logvar = torch.zeros(n_nodes, 3)
    preds = torch.zeros(4)
    labels = torch.zeros(4)
    loss = gcn_loss(preds, labels, mu, logvar, n_nodes=n_nodes, norm=1.0)
    # BCE(0 logits, 0 labels) = log(2) != 0, but KLD contribution should vanish;
    # verify by comparing to the BCE term alone.
    bce_only = torch.nn.functional.binary_cross_entropy_with_logits(preds, labels)
    assert torch.isclose(loss, bce_only, atol=1e-6)


def test_gcn_loss_increases_away_from_prior():
    n_nodes = 5
    mu_prior = torch.zeros(n_nodes, 3)
    logvar_prior = torch.zeros(n_nodes, 3)
    mu_shifted = torch.ones(n_nodes, 3) * 5
    preds = torch.zeros(4)
    labels = torch.zeros(4)

    loss_prior = gcn_loss(preds, labels, mu_prior, logvar_prior, n_nodes=n_nodes, norm=1.0)
    loss_shifted = gcn_loss(preds, labels, mu_shifted, logvar_prior, n_nodes=n_nodes, norm=1.0)
    assert loss_shifted > loss_prior
