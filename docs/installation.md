# Installation

NeuroFate targets Python 3.11 or 3.12 on Apple Silicon and Linux workstations.

```bash
conda env create -f environment.yml
conda activate neurofate
python -m pip install -e .
```

Optional MPS neural modeling requires PyTorch:

```bash
python -m pip install -e ".[torch]"
```

The package exposes a CLI after installation:

```bash
neurofate check-system
neurofate doctor
```
