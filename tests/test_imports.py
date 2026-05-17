from neurofate import __version__
from neurofate.features import planned_feature_groups
from neurofate.models import training_enabled


def test_version_is_defined():
    assert __version__


def test_planned_feature_groups_include_evolution():
    assert "evolutionary_conservation" in planned_feature_groups()


def test_training_disabled_in_skeleton():
    assert training_enabled() is False
