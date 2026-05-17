# Quickstart

Start with the lightweight checks:

```bash
neurofate check-system
neurofate doctor
python scripts/53_validate_neurofate_outputs.py
```

Generate an end-user report from existing outputs:

```bash
python scripts/51_generate_end_user_report.py
python scripts/52_generate_reproducibility_manifest.py
python scripts/54_no_overclaiming_audit.py
```

NeuroFate does not download datasets or run heavy analysis automatically.
