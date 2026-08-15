import pytest
import yaml

from graphist.config import DatasetConfig


MINIMAL_YAML = """
name: toy_dataset
data:
  loading_mode: manual
  counts_file: /tmp/counts.txt
  meta_file: /tmp/meta.txt
train:
  epochs: 5
"""


def test_from_yaml_loads_and_applies_defaults(tmp_path):
    path = tmp_path / "toy.yaml"
    path.write_text(MINIMAL_YAML)
    config = DatasetConfig.from_yaml(str(path))

    assert config.name == "toy_dataset"
    assert config.data.counts_file == "/tmp/counts.txt"
    assert config.train.epochs == 5
    # unset fields fall back to dataclass defaults
    assert config.train.lr == 1e-4
    assert config.model.k_neighbors == 12


def test_from_dict_requires_name():
    with pytest.raises(ValueError):
        DatasetConfig.from_dict({"data": {"loading_mode": "manual"}})


def test_real_configs_load_without_error():
    """The 4 shipped dataset configs must at least parse cleanly."""
    import os

    configs_dir = os.path.join(os.path.dirname(__file__), "..", "configs")
    for fname in ["pdac.yaml", "brca_pacsi.yaml", "brca_commot.yaml", "maynard.yaml"]:
        config = DatasetConfig.from_yaml(os.path.join(configs_dir, fname))
        assert config.name
        assert config.model.gmt_path
