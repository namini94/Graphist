"""Per-dataset configuration.

Every value here used to be a hardcoded constant or absolute
``/Users/naminiyakan/...`` path baked into one of the four original pipeline
scripts (see ``stage2_pathway_vgae/legacy/``). Defaults on each field match
the *PDAC* script (the fullest/most-recently-maintained one); every other
dataset's config file overrides only the fields that actually differ,
per the diff table in the refactor plan.
"""
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import yaml


@dataclass
class DataConfig:
    loading_mode: str = "manual"  # "manual" (raw counts + meta txt) or "visium" (10x Visium dir)
    data_root: Optional[str] = None  # visium mode: path to the 10x outs/ directory
    counts_file: Optional[str] = None  # manual mode: raw counts file (genes x spots)
    meta_file: Optional[str] = None  # manual mode: metadata file (spot coords + annotations)
    annotation_file: Optional[str] = None  # visium mode: optional extra annotation table to join
    annotation_sep: str = "\t"
    annotation_column: str = "annotations"  # obs column used for ARI / boxplot grouping
    annotation_column_source: Optional[str] = None  # if set, copy obs[source] -> obs[annotation_column]
    annotation_relabel: Dict[str, str] = field(default_factory=dict)  # raw value -> display value

    # Optional Stage-1 (Scissor-style regression) spot-selection CSV, joined in as an
    # extra obs column. Used by PDAC/BRCA-PACSI to group spots into
    # Graphist(+)/Graphist(-)/Background for differential pathway-activity testing.
    selection_file: Optional[str] = None
    selection_sep: str = ","
    selection_column: str = "selection"
    selection_relabel: Dict[int, str] = field(
        default_factory=lambda: {0: "Background", 1: "Graphist (+)", 2: "Graphist (-)"}
    )


@dataclass
class PreprocessConfig:
    min_cells: int = 20
    min_counts: int = 10
    target_sum: float = 1e6
    n_top_genes: int = 5000
    use_count_layer: bool = False
    n_pca_components: int = 200
    pca_seed: int = 42


@dataclass
class ModelConfig:
    gmt_path: str = ""
    add_nodes: int = 1
    min_genes: int = 0
    max_genes: int = 5000
    gcn_hidden1: int = 800
    p_drop: float = 0.2
    dropout: float = 0.1
    positive_decoder: bool = True
    k_neighbors: int = 12


@dataclass
class TrainConfig:
    seed: int = 2024
    epochs: int = 200
    lr: float = 1e-4
    decay: float = 5e-4
    rec_w: float = 1000
    gcn_w: float = 0.1
    gcn_rec_w: float = 10


@dataclass
class AnalysisConfig:
    output_dir: str = "."
    mclust_labels_path: Optional[str] = None
    pathways_to_plot: List[str] = field(default_factory=list)
    boxplot_pathways: List[str] = field(default_factory=list)
    boxplot_palette: Dict[str, str] = field(default_factory=dict)
    de_enabled: bool = False
    de_group_column: Optional[str] = None
    de_group1: Optional[str] = None
    de_group2: Optional[str] = None
    de_mode: str = "vanilla"


@dataclass
class DatasetConfig:
    name: str
    data: DataConfig = field(default_factory=DataConfig)
    preprocess: PreprocessConfig = field(default_factory=PreprocessConfig)
    model: ModelConfig = field(default_factory=ModelConfig)
    train: TrainConfig = field(default_factory=TrainConfig)
    analysis: AnalysisConfig = field(default_factory=AnalysisConfig)

    @classmethod
    def from_yaml(cls, path: str) -> "DatasetConfig":
        with open(path) as f:
            raw = yaml.safe_load(f)
        return cls.from_dict(raw)

    @classmethod
    def from_dict(cls, raw: dict) -> "DatasetConfig":
        if "name" not in raw:
            raise ValueError("Dataset config must have a 'name' field")
        return cls(
            name=raw["name"],
            data=DataConfig(**raw.get("data", {})),
            preprocess=PreprocessConfig(**raw.get("preprocess", {})),
            model=ModelConfig(**raw.get("model", {})),
            train=TrainConfig(**raw.get("train", {})),
            analysis=AnalysisConfig(**raw.get("analysis", {})),
        )
