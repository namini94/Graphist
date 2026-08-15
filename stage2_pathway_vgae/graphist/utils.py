"""General utilities shared across the graphist package."""
import os
import random

import numpy as np
import torch


def seed_everything(seed: int = 2024) -> None:
    """Seed python, numpy, and torch RNGs (plus cudnn determinism flags).

    All four original per-dataset scripts used this seed value in some form,
    but two of them (Maynard, BRCA-COMMOT) never actually called it, making
    those runs non-reproducible. The refactored pipeline always calls this.
    """
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
