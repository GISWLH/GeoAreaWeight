"""Write docs/summary_table.csv and docs/summary_table.md from docs/results.json."""
import csv
import json
from pathlib import Path

DOCS = Path(__file__).resolve().parents[1] / "docs"
rows = json.load(open(DOCS / "results.json"))["rows"]


def fmt(v, ref):
    """Round like the rest of the table: more decimals for small magnitudes."""
    a = abs(ref)
    return f"{v:.3f}" if a < 1 else f"{v:.2f}" if a < 100 else f"{v:.1f}"


cols = ["dataset / statistic", "kind", "unit", "period", "grid", "mean_unweighted", "mean_area_weighted",
        "mean_diff (unw-w)", "mean_rel_diff_%", "trend_unweighted (per decade)", "trend_area_weighted",
        "trend_diff", "trend_rel_diff_%", "note"]
with open(DOCS / "summary_table.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(cols)
    for r in rows:
        w.writerow([r["name"], r["kind"], r["unit"], r["period"], r["grid"], f"{r['mean_unw']:.4f}",
                    f"{r['mean_w']:.4f}", f"{r['mean_diff']:.4f}", f"{r['mean_rel_pct']:.2f}",
                    f"{r['trend_unw']:.4f}", f"{r['trend_w']:.4f}", f"{r['trend_diff']:.4f}",
                    f"{r['trend_rel_pct']:.1f}", r["note"]])

md = ["| Dataset / region (unit) | Period | Mean, unweighted | Mean, area-weighted | Difference | Rel. diff. "
      "| Trend per decade, unweighted | Trend, area-weighted | Rel. diff. |",
      "|---|---|---|---|---|---|---|---|---|"]
for r in rows:
    m, t = r["mean_w"], r["trend_w"]
    md.append(f"| {r['name']} ({r['unit']}) | {r['period']} | {fmt(r['mean_unw'], m)} | {fmt(m, m)} | "
              f"{fmt(r['mean_diff'], m)} | {r['mean_rel_pct']:.1f}% | {fmt(r['trend_unw'], t)} | {fmt(t, t)} | "
              f"{r['trend_rel_pct']:.0f}% |")
(DOCS / "summary_table.md").write_text("\n".join(md) + "\n", encoding="utf-8")
print("wrote", DOCS / "summary_table.csv", DOCS / "summary_table.md")
