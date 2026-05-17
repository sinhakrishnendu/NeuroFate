# Runtime Benchmarks

These values are current project observations or placeholders and should be refreshed for each release.

| Workflow | Observed or Expected Runtime |
| --- | --- |
| SEA-AD 34 GB download | User-dependent |
| Target-gene sparse extraction | About 10 minutes on M5 Max |
| Phase 3 sparse summaries | Under 1 minute |
| Phase 4 statistics | About 2 minutes |
| Phase 5 CPU ML | Under 1 minute |
| Phase 6 MPS MLP | Seconds on current donor table |
| Mathys CSV extraction | Seconds |

Runtime depends on storage speed, available memory, Python environment, and exact cohort files. NeuroFate avoids loading full H5AD expression matrices into memory.
