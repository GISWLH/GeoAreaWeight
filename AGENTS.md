# Notes for contributors and coding agents

Users of the library should read `README.md` (or `llms.txt` for a compact API
summary). This file is about changing the code.

## Layout

- `python/src/geoareaweight/weights.py`: weights, bounds, cell areas (NumPy only).
- `python/src/geoareaweight/means.py`: `area_mean`, `region_mean`, `area_integral`,
  `area_fraction`; NumPy and xarray dispatch (`_dispatch`), reduction (`_reduce_core`).
- `R/R/geoareaweight.R`, `R/man/geoareaweight.Rd`: R package (hand-written Rd).
- `matlab/gaw_*.m` public functions, `matlab/private/` internal helpers.
- `validation/`: the same 312 numbers computed in all three languages.

## Rules

- **Keep the three implementations in sync.** A behaviour change in one language
  needs the same change (and test) in the others, plus regenerated
  `validation/ref_*.csv` and a passing `validation/compare.py`.
- Angles are degrees; 2-D fields are `(nlat, nlon)`; default method `"band"`.
- Do not change numerical results silently: `validation/ref_python.csv` is the
  reference; differences above ~1e-13 relative need an explanation in `CHANGELOG.md`.
- Error messages should say what to do (e.g. "pass lon_bounds").
- Public Python functions need type hints and numpydoc docstrings; examples in
  docstrings are run as doctests.

## Commands

```bash
cd python && pip install -e ".[dev]" && ruff check . && pytest      # Python
R CMD build R && R CMD check --no-manual geoareaweight_*.tar.gz    # R (repo root)
octave-cli --eval "addpath('matlab'); addpath('matlab/tests'); test_gaw"
cd validation && python ref_python.py ref_python.csv && Rscript ref_r.R ref_r.csv \
  && octave-cli --eval "addpath('../matlab'); ref_octave('ref_octave.csv')" && python compare.py
```

## Releasing

Bump the version in `python/src/geoareaweight/__init__.py`, `R/DESCRIPTION`,
`matlab/Contents.m` and `CITATION.cff`, update `CHANGELOG.md`, tag `vX.Y.Z` and
publish a GitHub release; `.github/workflows/release.yml` uploads the Python
package to PyPI (trusted publishing must be configured on PyPI first).
