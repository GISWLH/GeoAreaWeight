"""Replace the Chinese labels of the original fig4 rendering with English ones.

Used once (v0.1.0) because the time series behind fig4 were not stored in the
repository. Once docs/series.json exists, make_figs.py regenerates fig4 from
data and this script is no longer needed. The plotted data are not modified:
only text areas outside the data lines are blanked and re-labelled.

    python experiments/relabel_fig4.py ORIGINAL_PNG OUTPUT_PNG
"""
import sys

import matplotlib.font_manager as fm
import numpy as np
from PIL import Image, ImageDraw, ImageFont

src, dst = sys.argv[1], sys.argv[2]
img = Image.open(src).convert("RGB")
S = img.width / 2000.0                       # boxes below are in 2000-px-wide display coordinates
font_path = fm.findfont("DejaVu Sans")
arr = np.asarray(img)


def ink_box(box, thr=200):
    """Tight bounding box of non-white pixels inside a display-coordinate box."""
    x0, y0, x1, y1 = (int(round(v * S)) for v in box)
    sub = arr[y0:y1, x0:x1].min(axis=2) < thr
    ys, xs = np.nonzero(sub)
    return x0 + xs.min(), y0 + ys.min(), x0 + xs.max() + 1, y0 + ys.max() + 1


def font(px):
    return ImageFont.truetype(font_path, int(px))


d = ImageDraw.Draw(img)


def replace_centered(box, lines, size=38, pad=4):
    x0, y0, x1, y1 = ink_box(box)
    d.rectangle([x0 - pad, y0 - pad, x1 + pad, y1 + pad], fill="white")
    cx = (x0 + x1) / 2
    lh = (y1 - y0) / len(lines)
    for i, t in enumerate(lines):
        d.text((cx, y0 + lh * (i + 0.5)), t, font=font(size), fill="black", anchor="mm")


def replace_left(box, text, size=38, pad=4):
    x0, y0, x1, y1 = ink_box(box)
    d.rectangle([x0 - pad, y0 - pad, x1 + pad, y1 + pad], fill="white")
    # white halo keeps the text legible where a data line passes behind it
    d.text((x0, (y0 + y1) / 2), text, font=font(size), fill="black", anchor="lm",
           stroke_width=6, stroke_fill="white")


def replace_vertical(box, text, size=38, pad=4):
    x0, y0, x1, y1 = ink_box(box)
    d.rectangle([x0 - pad, y0 - pad, x1 + pad, y1 + pad], fill="white")
    f = font(size)
    w = int(d.textlength(text, font=f)) + 8
    tmp = Image.new("RGB", (w, size + 16), "white")
    ImageDraw.Draw(tmp).text((4, (size + 16) / 2), text, font=f, fill="black", anchor="lm")
    tmp = tmp.rotate(90, expand=True)
    img.paste(tmp, (int((x0 + x1) / 2 - tmp.width / 2), int((y0 + y1) / 2 - tmp.height / 2)))


# panel titles
replace_centered((190, 8, 900, 70), [
    "(a) NCEP/NCAR R1 near-surface air temperature, global anomaly (vs 1981–2010)",
    "trend: unweighted 0.33 vs weighted 0.19 °C/decade (+75%)"])
replace_centered((1200, 8, 1880, 70), [
    "(b) GISTEMP v4 (2°) global temperature anomaly, 1880–2020",
    "1980–2020 trend: unweighted 0.26 vs weighted 0.20 °C/decade (+34%)"])
replace_centered((180, 660, 925, 722), [
    "(c) Land area fraction with SPEI-12 ≤ −1 (CSIC SPEIbase v2.10, 0.5°)",
    "mean 21.4% vs 21.9%; trend 2.80 vs 2.83 pp/decade: small difference"])
replace_centered((1190, 660, 1900, 722), [
    "(d) GPCC v2020 annual land precipitation (1°)",
    "mean 774 vs 902 mm/yr (−14%); trend 5.6 vs 5.6 mm/yr per decade"])
# legends (keep the line handles, replace the text)
for x, y in [(200, 95), (1195, 95), (200, 745), (1195, 745)]:
    replace_left((x, y, x + 110, y + 25), "Unweighted", size=40)
    replace_left((x, y + 35, x + 110, y + 62), "Area-weighted", size=40)
# axis labels
replace_vertical((12, 270, 48, 375), "Anomaly (°C)", size=40)
replace_vertical((1008, 175, 1042, 465), "Anomaly (°C, vs 1981–2010)", size=40)
replace_vertical((18, 900, 52, 1050), "Area fraction (%)", size=40)
for box in [(530, 615, 565, 645), (1530, 615, 1565, 645), (535, 1265, 570, 1297), (1530, 1265, 1565, 1297)]:
    replace_centered(box, ["Year"], size=40)
img.save(dst, optimize=True)
print("wrote", dst)
