from pathlib import Path
import importlib.util
import sys
import types


SCRIPT = Path("scripts/13_decode_sea_ad_categories.py")


def test_category_decoder_avoids_scanpy_and_matrix_loaders():
    text = SCRIPT.read_text(encoding="utf-8")
    lowered = text.lower()
    assert "import scanpy" not in lowered
    assert "read_h5ad" not in lowered
    assert "anndata" not in lowered
    assert "AnnData" not in text


def test_category_decoder_has_x_access_guard_without_direct_x_indexing():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "FORBIDDEN_ROOT_KEY" in text
    assert "RuntimeError" in text
    assert '["X"]' not in text
    assert "['X']" not in text
    assert "obs/__categories" in text


def test_category_decoder_declares_expected_outputs():
    text = SCRIPT.read_text(encoding="utf-8")
    for filename in [
        "sea_ad_obs_metadata_decoded.tsv",
        "table1_sea_ad_publication_ready.tsv",
        "sea_ad_category_mapping.tsv",
        "sea_ad_unsupported_category_nodes.tsv",
        "13_decode_sea_ad_categories.log",
    ]:
        assert filename in text


def test_category_decoder_has_no_analysis_workflows():
    text = SCRIPT.read_text(encoding="utf-8")
    for forbidden in ["normalize_total", "pca", "umap", "leiden", "neighbors", "fit("]:
        assert forbidden not in text


def load_decoder_module(monkeypatch=None):
    if monkeypatch is not None:
        monkeypatch.setitem(
            sys.modules,
            "h5py",
            types.SimpleNamespace(Dataset=FakeDataset, Group=FakeGroup, File=object),
        )
    spec = importlib.util.spec_from_file_location("decode_sea_ad_categories", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class FakeArray(list):
    @property
    def shape(self):
        return (len(self),)


class FakeDataset:
    def __init__(self, name, values):
        self.name = name
        self._values = FakeArray(values)

    def __getitem__(self, key):
        assert key == ()
        return self._values


class FakeGroup(dict):
    def __init__(self, name):
        super().__init__()
        self.name = name


def test_unsupported_category_nodes_are_nonfatal(monkeypatch):
    module = load_decoder_module(monkeypatch)
    root = FakeGroup("/")
    obs = FakeGroup("/obs")
    categories = FakeGroup("/obs/__categories")
    categories["Class"] = FakeDataset("/obs/__categories/Class", [b"Neuronal", b"Glial"])
    categories["Hispanic"] = FakeGroup("/obs/__categories/Hispanic")
    obs["__categories"] = categories
    root["obs"] = obs

    reader = module.SafeH5Reader(root)
    mappings, unsupported = module.load_category_mappings(reader)

    assert mappings["Class"] == {"0": "Neuronal", "1": "Glial"}
    assert len(unsupported) == 1
    assert unsupported[0]["category_column"] == "Hispanic"
    assert unsupported[0]["node_type"] == "Group"
