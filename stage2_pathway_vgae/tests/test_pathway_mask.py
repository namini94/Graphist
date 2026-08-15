import numpy as np

from graphist.data.pathway_mask import create_mask, pathway_names, read_gmt, write_gmt


def test_read_gmt_round_trip(tmp_path):
    original = {"PATH_1": ["G1", "G2", "G3"], "PATH_2": ["G4", "G5"]}
    path = tmp_path / "roundtrip.gmt"
    write_gmt(original, str(path))
    loaded = read_gmt(str(path))
    assert dict(loaded) == original


def test_read_gmt_size_filters(tmp_path):
    original = {"SMALL": ["G1"], "BIG": [f"G{i}" for i in range(10)]}
    path = tmp_path / "filtered.gmt"
    write_gmt(original, str(path))
    loaded = read_gmt(str(path), min_g=2, max_g=20)
    assert "SMALL" not in loaded
    assert "BIG" in loaded


def test_create_mask_shape_and_unannotated_columns(tiny_gmt_path):
    features = [f"GENE{i}" for i in range(60)]
    mask = create_mask(features, tiny_gmt_path, add_nodes=2)
    # 5 toy pathways + 2 unannotated nodes
    assert mask.shape == (60, 7)
    # last 2 columns are the fully-connected unannotated nodes
    np.testing.assert_array_equal(mask[:, -2:], np.ones((60, 2)))


def test_create_mask_membership_correctness(tiny_gmt_path):
    features = [f"GENE{i}" for i in range(60)]
    mask = create_mask(features, tiny_gmt_path, add_nodes=1)
    # PATHWAY_A = GENE0..GENE9 -> column 0
    assert mask[0, 0] == 1  # GENE0 in PATHWAY_A
    assert mask[15, 0] == 0  # GENE15 not in PATHWAY_A
    # PATHWAY_B = GENE5..GENE14 overlaps PATHWAY_A on GENE5-9
    assert mask[7, 0] == 1 and mask[7, 1] == 1


def test_pathway_names_matches_mask_columns(tiny_gmt_path):
    features = [f"GENE{i}" for i in range(60)]
    mask = create_mask(features, tiny_gmt_path, add_nodes=1)
    names = pathway_names(tiny_gmt_path, add_nodes=1)
    assert len(names) == mask.shape[1]
    assert names[-1] == "UNANNOTATED_0"
    assert names[:5] == ["PATHWAY_A", "PATHWAY_B", "PATHWAY_C", "PATHWAY_D", "PATHWAY_E"]
