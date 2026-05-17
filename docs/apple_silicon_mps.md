# Apple Silicon MPS

NeuroFate includes a small donor-level PyTorch MPS MLP for Apple Silicon.

Check device support:

```bash
python scripts/26_check_mps_device.py
```

The MPS model is intentionally small and donor-level. It is not a transformer and is not a large model.
