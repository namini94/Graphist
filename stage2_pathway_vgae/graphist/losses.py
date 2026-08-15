"""Loss functions for the pathway-activity VGAE."""
import torch
import torch.nn.functional as F


def reconstruction_loss(decoded: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
    """Mean-squared-error reconstruction loss."""
    return torch.nn.MSELoss()(decoded, x)


def gcn_loss(
    preds: torch.Tensor,
    labels: torch.Tensor,
    mu: torch.Tensor,
    logvar: torch.Tensor,
    n_nodes: int,
    norm: float,
) -> torch.Tensor:
    """VGAE loss: weighted edge-reconstruction BCE + analytic Gaussian KL divergence.

    See Appendix B of Kingma & Welling, "Auto-Encoding Variational Bayes", ICLR 2014
    (https://arxiv.org/abs/1312.6114): KL = -0.5 * sum(1 + log(sigma^2) - mu^2 - sigma^2).
    Note ``logvar`` here is actually the encoder's raw second output (not literally
    log-variance) — squared before exponentiating, matching the original code exactly.
    """
    cost = norm * F.binary_cross_entropy_with_logits(preds, labels)
    kld = -0.5 / n_nodes * torch.mean(
        torch.sum(1 + 2 * logvar - mu.pow(2) - logvar.exp().pow(2), 1)
    )
    return cost + kld
