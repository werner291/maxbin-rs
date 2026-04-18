# Versioning

maxbin-rs uses version numbers inspired by [Semantic Versioning](https://semver.org/),
but adapted to a bioinformatics tool where the API contract is the binning
output, not a programming interface. **"Breaking change" means "changes
binning output on the same input."**

A crash fix, a better error message, or a new CLI flag is not breaking.
A different default threshold that changes which contigs end up in which
bin is.

## Version series

### v0.1.x — Bug-for-bug compatible

Output matches the original MaxBin2 (modulo seed ordering, which is
non-deterministic in both tools). This is the "drop-in replacement" series.

Allowed changes:
- Crash fixes, error handling, better diagnostics
- Packaging and build fixes
- Performance improvements that don't change output
- Documentation

Not allowed:
- Anything that changes bin assignments on the same input

### v0.2.x — Corrected defaults, obvious bug fixes

Same algorithm, but with corrected defaults and fixed obvious bugs.
A user can recover v0.1 behavior by passing the original's flags
explicitly (e.g. `-prob_threshold 0.5`).

Changes so far:
- `prob_threshold` default: 0.5 → 0.9 (matching the documented value)

Each breaking change is documented with motivation — why the old
behavior was wrong, not just what changed.

### v0.3+ — New features

Prodigal gene caller, GTDB markers, parallelized outer loop, etc.
Version numbers TBD when we get there.

## Branching

`main` should stay in a working state. Work-in-progress goes on feature
branches and merges to main when it's ready (`--ff-only`). No PR
ceremony required — the merge is the review.

When main moves past a version series (e.g., starts v0.2 work), a
`release/0.1` branch is created from the last v0.1.x tag for any
necessary backports. Tags mark releases.

There is currently one maintainer. If contributors appear, this scales
to PR-based review naturally — the versioning contract stays the same
regardless of process.

