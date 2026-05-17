#!/usr/bin/env python3
"""Check PyTorch Metal/MPS availability for NeuroFate donor-level modeling."""

from __future__ import annotations


def main() -> int:
    try:
        import torch
    except ImportError:
        print("torch version: not installed")
        print("MPS built: false")
        print("MPS available: false")
        print("selected device: cpu")
        print("small tensor test result: skipped because torch is not installed")
        return 1

    mps_built = bool(torch.backends.mps.is_built())
    mps_available = bool(torch.backends.mps.is_available())
    device = "mps" if mps_available else "cpu"
    tensor = torch.tensor([1.0, 2.0, 3.0], device=device)
    result = float((tensor * 2.0).sum().detach().cpu().item())

    print(f"torch version: {torch.__version__}")
    print(f"MPS built: {mps_built}")
    print(f"MPS available: {mps_available}")
    print(f"selected device: {device}")
    print(f"small tensor test result: {result}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
