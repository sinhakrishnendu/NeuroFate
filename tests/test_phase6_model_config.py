from pathlib import Path


CONFIG = Path("configs/neurofate_mps_model_config.yaml")


def parse_simple_yaml(path: Path) -> dict[str, object]:
    parsed: dict[str, object] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        key, value = line.split(":", 1)
        value = value.strip()
        if value.replace(".", "", 1).isdigit():
            parsed[key.strip()] = float(value) if "." in value else int(value)
        else:
            parsed[key.strip()] = value
    return parsed


def test_phase6_mps_config_exists_and_has_required_keys():
    config = parse_simple_yaml(CONFIG)
    expected = {
        "seed",
        "device_preference",
        "epochs",
        "batch_size",
        "hidden_dim",
        "dropout",
        "learning_rate",
        "weight_decay",
        "early_stopping_patience",
        "test_size",
        "validation_size",
    }
    assert expected.issubset(config)


def test_phase6_mps_config_values_match_requested_defaults():
    config = parse_simple_yaml(CONFIG)
    assert config["seed"] == 42
    assert config["device_preference"] == "mps"
    assert config["epochs"] == 300
    assert config["batch_size"] == 16
    assert config["hidden_dim"] == 64
    assert config["dropout"] == 0.25
    assert config["learning_rate"] == 0.001
    assert config["weight_decay"] == 0.0001
    assert config["early_stopping_patience"] == 40
    assert config["test_size"] == 0.25
    assert config["validation_size"] == 0.25
