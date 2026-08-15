"""graphist: a spatially-aware, pathway-masked variational graph autoencoder
for quantifying gene-pathway activity in spatial transcriptomics data.

See ``graphist.pipeline.run_pipeline`` for the high-level entrypoint, or use
the individual pieces (``graphist.model``, ``graphist.trainer``,
``graphist.data``) directly for custom workflows.
"""
from .config import DatasetConfig
from .model.graphist_model import GraphistModel
from .trainer import GraphistTrainer

__all__ = ["DatasetConfig", "GraphistModel", "GraphistTrainer"]
__version__ = "0.1.0"
