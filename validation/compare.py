"""Compare the reference values computed by the Python, R and Octave implementations.

    python compare.py                                   # ref_python.csv vs ref_r.csv, ref_octave.csv
    python compare.py PY.csv R.csv OCTAVE.csv [COMMITTED_PY.csv]

Each CSV holds lines ``key[i],value``. Exits with status 1 if any pair differs by
more than a relative tolerance of 1e-12 (NaN must match NaN).
"""
import math
import sys

TOL = 1e-12


def load(path):
    out = {}
    with open(path) as f:
        for line in f:
            if line.strip():
                k, v = line.strip().rsplit(",", 1)
                out[k] = float("nan" if v in ("NA", "NaN", "nan") else v)
    return out


def compare(ref, other, label):
    if ref.keys() != other.keys():
        print(f"{label}: different keys: {sorted(set(ref) ^ set(other))[:10]}")
        return False
    worst, where = 0.0, None
    for k, a in ref.items():
        b = other[k]
        if math.isnan(a) or math.isnan(b):
            if not (math.isnan(a) and math.isnan(b)):
                print(f"{label}: NaN mismatch at {k}")
                return False
            continue
        r = abs(a - b) / max(abs(a), 1e-300) if a != 0 else abs(b)
        if r > worst:
            worst, where = r, k
    print(f"{label}: {len(ref)} numbers, max rel diff {worst:.2e} at {where}")
    return worst <= TOL


args = sys.argv[1:] or ["ref_python.csv", "ref_r.csv", "ref_octave.csv"]
py = load(args[0])
ok = compare(py, load(args[1]), "python vs R") & compare(py, load(args[2]), "python vs octave")
if len(args) > 3:
    ok &= compare(load(args[3]), py, "committed reference vs python")
sys.exit(0 if ok else 1)
