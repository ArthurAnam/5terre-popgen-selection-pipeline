#!/usr/bin/env python3
import argparse
import csv
import hashlib
import html
import math
import os
from pathlib import Path

POP_ORDER = ["CT", "CEU", "FIN", "GBR", "IBS", "TSI"]
POP_COLORS = {
    "CT": "#000000",
    "CEU": "#0072B2",
    "FIN": "#E69F00",
    "GBR": "#009E73",
    "IBS": "#D55E00",
    "TSI": "#CC79A7",
}
POP_SHAPES = {
    "CT": "star",
    "CEU": "circle",
    "FIN": "square",
    "GBR": "triangle",
    "IBS": "diamond",
    "TSI": "cross",
}

def read_tsv(path):
    with open(path, "r", encoding="utf-8") as h:
        return list(csv.DictReader(h, delimiter="\t"))

def read_registry(path):
    return {r["figure_id"]: r for r in read_tsv(path)}

def finite(values):
    return [float(v) for v in values if v not in ("", "NA", "nan", "NaN") and math.isfinite(float(v))]

def esc(x):
    return html.escape(str(x), quote=True)

def jitter(key, width=0.28):
    digest = hashlib.md5(str(key).encode("utf-8")).hexdigest()
    unit = int(digest[:8], 16) / 0xFFFFFFFF
    return (unit - 0.5) * 2.0 * width

def q(values, frac):
    x = sorted(values)
    if not x:
        return float("nan")
    if len(x) == 1:
        return x[0]
    pos = (len(x) - 1) * frac
    lo, hi = int(math.floor(pos)), int(math.ceil(pos))
    if lo == hi:
        return x[lo]
    w = pos - lo
    return x[lo] * (1 - w) + x[hi] * w

def bounds(values):
    vals = finite(values)
    if not vals:
        return 0.0, 1.0
    lo, hi = min(vals), max(vals)
    if lo == hi:
        pad = abs(lo) * 0.1 or 1.0
    else:
        pad = (hi - lo) * 0.07
    return lo - pad, hi + pad

def ticks(lo, hi, n=6):
    return [lo + i * (hi - lo) / (n - 1) for i in range(n)]

def fmt(v):
    av = abs(v)
    if av >= 1000 or (av > 0 and av < 0.001):
        return f"{v:.3e}"
    if av >= 10:
        return f"{v:.1f}"
    return f"{v:.4f}".rstrip("0").rstrip(".")

def axes_svg(width, height, ml, mr, mt, mb, xlo, xhi, ylo, yhi, xlabel, ylabel):
    pw, ph = width - ml - mr, height - mt - mb
    def sx(x):
        return ml + (x - xlo) / (xhi - xlo) * pw
    def sy(y):
        return mt + (yhi - y) / (yhi - ylo) * ph
    parts = [
        f'<line x1="{ml}" y1="{height-mb}" x2="{width-mr}" y2="{height-mb}" class="axis"/>',
        f'<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{height-mb}" class="axis"/>',
    ]
    for x in ticks(xlo, xhi):
        px = sx(x)
        parts.append(f'<line x1="{px:.2f}" y1="{height-mb}" x2="{px:.2f}" y2="{height-mb+6}" class="tick"/>')
        parts.append(f'<text x="{px:.2f}" y="{height-mb+22}" text-anchor="middle" class="ticktext">{esc(fmt(x))}</text>')
    for y in ticks(ylo, yhi):
        py = sy(y)
        parts.append(f'<line x1="{ml-6}" y1="{py:.2f}" x2="{ml}" y2="{py:.2f}" class="tick"/>')
        parts.append(f'<line x1="{ml}" y1="{py:.2f}" x2="{width-mr}" y2="{py:.2f}" class="grid"/>')
        parts.append(f'<text x="{ml-10}" y="{py+4:.2f}" text-anchor="end" class="ticktext">{esc(fmt(y))}</text>')
    parts.append(f'<text x="{ml+pw/2:.2f}" y="{height-18}" text-anchor="middle" class="label">{esc(xlabel)}</text>')
    parts.append(f'<text x="20" y="{mt+ph/2:.2f}" text-anchor="middle" class="label" transform="rotate(-90 20 {mt+ph/2:.2f})">{esc(ylabel)}</text>')
    return parts, sx, sy

def point_svg(x, y, pop, tip, size=5.0, opacity=0.80):
    color = POP_COLORS.get(pop, "#3366AA")
    shape = POP_SHAPES.get(pop, "circle")
    attr = f'class="hoverpoint" data-tooltip="{esc(tip)}"'
    if shape == "square":
        return f'<rect x="{x-size:.2f}" y="{y-size:.2f}" width="{2*size:.2f}" height="{2*size:.2f}" fill="{color}" fill-opacity="{opacity}" {attr}/>'
    if shape == "triangle":
        pts = f"{x:.2f},{y-size*1.15:.2f} {x-size*1.05:.2f},{y+size*.9:.2f} {x+size*1.05:.2f},{y+size*.9:.2f}"
        return f'<polygon points="{pts}" fill="{color}" fill-opacity="{opacity}" {attr}/>'
    if shape == "diamond":
        pts = f"{x:.2f},{y-size*1.2:.2f} {x-size:.2f},{y:.2f} {x:.2f},{y+size*1.2:.2f} {x+size:.2f},{y:.2f}"
        return f'<polygon points="{pts}" fill="{color}" fill-opacity="{opacity}" {attr}/>'
    if shape == "cross":
        return (f'<g {attr} stroke="{color}" stroke-width="2.2" stroke-linecap="round" opacity="{opacity}">'
                f'<line x1="{x-size:.2f}" y1="{y-size:.2f}" x2="{x+size:.2f}" y2="{y+size:.2f}"/>'
                f'<line x1="{x-size:.2f}" y1="{y+size:.2f}" x2="{x+size:.2f}" y2="{y-size:.2f}"/></g>')
    if shape == "star":
        pts = []
        for i in range(10):
            angle = -math.pi / 2 + i * math.pi / 5
            radius = size * (1.35 if i % 2 == 0 else 0.58)
            pts.append(f"{x + radius*math.cos(angle):.2f},{y + radius*math.sin(angle):.2f}")
        return f'<polygon points="{" ".join(pts)}" fill="{color}" fill-opacity="{opacity}" stroke="#000" stroke-width="0.7" {attr}/>'
    return f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{size:.2f}" fill="{color}" fill-opacity="{opacity}" {attr}/>'

def page(title, caption, svg, notes=""):
    return f"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{esc(title)}</title>
<style>
body{{font-family:Arial,Helvetica,sans-serif;margin:24px;color:#222;max-width:1180px}}
h1{{font-size:24px;margin-bottom:8px}} .caption{{font-size:15px;line-height:1.45;max-width:1000px}}
.note{{font-size:13px;color:#555}} svg{{width:100%;height:auto;border:1px solid #ddd;background:white}}
.axis,.tick{{stroke:#333;stroke-width:1}} .grid{{stroke:#ddd;stroke-width:0.7}} .ticktext{{font-size:11px;fill:#333}}
.label{{font-size:13px;fill:#222}} .legend{{font-size:12px;fill:#222}} .hoverpoint{{cursor:crosshair}}
#tooltip{{display:none;position:fixed;z-index:9999;pointer-events:none;background:rgba(20,20,20,.94);color:white;
padding:8px 10px;border-radius:5px;font-size:13px;line-height:1.35;max-width:430px;box-shadow:0 2px 8px rgba(0,0,0,.28)}}
</style></head><body>
<h1>{esc(title)}</h1><p class="caption">{esc(caption)}</p>
<div id="tooltip"></div>
{svg}
<p class="note">Passa il mouse sui punti/bin: il tooltip mostra campione, popolazione e valori sottostanti. Gli HTML servono per ispezione; PNG/PDF e tabelle restano gli output per il manoscritto.</p>
{notes}
<script>
const tt=document.getElementById('tooltip');
document.querySelectorAll('.hoverpoint').forEach(el=>{{
  el.addEventListener('mouseenter',()=>{{tt.textContent=el.dataset.tooltip||'';tt.style.display='block';}});
  el.addEventListener('mousemove',(e)=>{{tt.style.left=(e.clientX+14)+'px';tt.style.top=(e.clientY+14)+'px';}});
  el.addEventListener('mouseleave',()=>{{tt.style.display='none';}});
}});
</script>
</body></html>"""

def scatter_html(rows, xfield, yfield, xlabel, ylabel, title, caption, tooltip_fields, popfield=None):
    width, height, ml, mr, mt, mb = 920, 610, 82, 40, 25, 62
    xs = [float(r[xfield]) for r in rows]
    ys = [float(r[yfield]) for r in rows]
    xlo, xhi = bounds(xs)
    ylo, yhi = bounds(ys)
    parts, sx, sy = axes_svg(width, height, ml, mr, mt, mb, xlo, xhi, ylo, yhi, xlabel, ylabel)
    for r in rows:
        pop = r.get(popfield, "") if popfield else ""
        tip = " | ".join(f"{k}: {r.get(k,'')}" for k in tooltip_fields)
        radius = 6.0 if pop == "CT" else 4.6
        parts.append(point_svg(sx(float(r[xfield])), sy(float(r[yfield])), pop, tip, radius, 0.82))
    if popfield:
        lx = width - mr - 120
        ly = mt + 12
        shown = [p for p in POP_ORDER if any(r.get(popfield) == p for r in rows)]
        for i, pop in enumerate(shown):
            y = ly + i * 22
            parts.append(point_svg(lx, y, pop, f"population: {pop}", 5.5, 1.0))
            parts.append(f'<text x="{lx+14}" y="{y+4}" class="legend">{pop}</text>')
    svg = f'<svg viewBox="0 0 {width} {height}" role="img">{"".join(parts)}</svg>'
    return page(title, caption, svg)

def grouped_html(rows, valuefield, ylabel, title, caption):
    width, height, ml, mr, mt, mb = 920, 610, 82, 35, 25, 62
    vals = [float(r[valuefield]) for r in rows]
    ylo, yhi = bounds(vals)
    ylo = min(0.0, ylo)
    pw, ph = width - ml - mr, height - mt - mb
    def sy(y):
        return mt + (yhi - y) / (yhi - ylo) * ph
    parts = [
        f'<line x1="{ml}" y1="{height-mb}" x2="{width-mr}" y2="{height-mb}" class="axis"/>',
        f'<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{height-mb}" class="axis"/>',
    ]
    for y in ticks(ylo, yhi):
        py = sy(y)
        parts.append(f'<line x1="{ml}" y1="{py:.2f}" x2="{width-mr}" y2="{py:.2f}" class="grid"/>')
        parts.append(f'<text x="{ml-10}" y="{py+4:.2f}" text-anchor="end" class="ticktext">{esc(fmt(y))}</text>')
    cats = [p for p in POP_ORDER if any(r["population"] == p for r in rows)]
    step = pw / len(cats)
    for i, pop in enumerate(cats):
        cx = ml + (i + 0.5) * step
        sub = [r for r in rows if r["population"] == pop]
        x = [float(r[valuefield]) for r in sub]
        vmin, vq1, vmed, vq3, vmax = min(x), q(x, .25), q(x, .5), q(x, .75), max(x)
        boxw = min(54, step * 0.45)
        color = POP_COLORS[pop]
        parts.append(f'<line x1="{cx}" y1="{sy(vmin):.2f}" x2="{cx}" y2="{sy(vmax):.2f}" stroke="{color}" stroke-width="1.4"/>')
        parts.append(f'<rect x="{cx-boxw/2:.2f}" y="{sy(vq3):.2f}" width="{boxw:.2f}" height="{sy(vq1)-sy(vq3):.2f}" fill="{color}" fill-opacity="0.13" stroke="{color}" stroke-width="1.5"/>')
        parts.append(f'<line x1="{cx-boxw/2:.2f}" y1="{sy(vmed):.2f}" x2="{cx+boxw/2:.2f}" y2="{sy(vmed):.2f}" stroke="{color}" stroke-width="2.4"/>')
        for r in sub:
            px = cx + jitter(r["IID"], min(24, step * 0.25))
            tip = (
                f"IID: {r['IID']} | population: {pop} | {valuefield}: {r[valuefield]} | "
                f"N_ROH>=1.5Mb: {r.get('N_ROH_GE1_5MB','')} | "
                f"total_ROH_Mb>=1.5: {r.get('TOTAL_ROH_MB_GE1_5MB','')} | "
                f"max_ROH_Mb: {r.get('MAX_ROH_MB','')}"
            )
            parts.append(point_svg(px, sy(float(r[valuefield])), pop, tip, 4.8 if pop != "CT" else 5.8, 0.76))
        parts.append(f'<text x="{cx:.2f}" y="{height-mb+23}" text-anchor="middle" class="ticktext">{pop}</text>')
    parts.append(f'<text x="20" y="{mt+ph/2:.2f}" text-anchor="middle" class="label" transform="rotate(-90 20 {mt+ph/2:.2f})">{esc(ylabel)}</text>')
    svg = f'<svg viewBox="0 0 {width} {height}" role="img">{"".join(parts)}</svg>'
    return page(title, caption, svg)

def ld_html(bins, summary, title, caption):
    rows = [r for r in bins if int(r["n_pairs"]) > 0 and r["mean_r2"] not in ("nan", "NaN", "")]
    width, height, ml, mr, mt, mb = 940, 600, 82, 35, 25, 62
    xs = [float(r["bin_center_kb"]) for r in rows]
    ys = [float(r["mean_r2"]) for r in rows]
    xlo, xhi = 0.0, max(xs)
    ylo, yhi = 0.0, max(ys) * 1.08
    parts, sx, sy = axes_svg(width, height, ml, mr, mt, mb, xlo, xhi, ylo, yhi, "Physical distance between SNPs (kb)", "Mean pairwise r²")
    pts = " ".join(f"{sx(float(r['bin_center_kb'])):.2f},{sy(float(r['mean_r2'])):.2f}" for r in rows)
    parts.append(f'<polyline points="{pts}" fill="none" stroke="#2F4B7C" stroke-width="1.8"/>')
    sm = {r["metric"]: r["value"] for r in summary}
    threshold = float(sm["decay_threshold_r2"])
    crossing = float(sm["first_crossing_bin_center_kb"])
    parts.append(f'<line x1="{ml}" y1="{sy(threshold):.2f}" x2="{width-mr}" y2="{sy(threshold):.2f}" stroke="#B22222" stroke-dasharray="7 5"/>')
    parts.append(f'<line x1="{sx(crossing):.2f}" y1="{mt}" x2="{sx(crossing):.2f}" y2="{height-mb}" stroke="#555" stroke-dasharray="3 5"/>')
    for r in rows:
        tip = f"distance bin center: {r['bin_center_kb']} kb | mean r2: {r['mean_r2']} | n_pairs: {r['n_pairs']} | SE: {r['se_r2']}"
        parts.append(f'<circle class="hoverpoint" data-tooltip="{esc(tip)}" cx="{sx(float(r["bin_center_kb"])):.2f}" cy="{sy(float(r["mean_r2"])):.2f}" r="3.0" fill="#2F4B7C" fill-opacity="0.85"/>')
    svg = f'<svg viewBox="0 0 {width} {height}" role="img">{"".join(parts)}</svg>'
    notes = (
        f'<p class="note">Baseline mean r² = {esc(sm["baseline_mean_r2"])}; '
        f'one-third threshold = {esc(sm["decay_threshold_r2"])}; '
        f'first crossing = {esc(sm["first_crossing_bin_center_kb"])} kb; '
        f'stable crossing = {esc(sm["stable_crossing_bin_center_kb"])} kb.</p>'
    )
    return page(title, caption, svg, notes)

def read_eigen(path):
    rows = read_tsv(path)
    return {int(r["PC"]): float(r["variance_percent"]) for r in rows}

def write(path, text):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(text, encoding="utf-8")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--registry", required=True)
    p.add_argument("--king-pairs", required=True)
    p.add_argument("--pca-unmasked", required=True)
    p.add_argument("--pca-unmasked-eigen", required=True)
    p.add_argument("--pca-masked", required=True)
    p.add_argument("--pca-masked-eigen", required=True)
    p.add_argument("--roh-individual", required=True)
    p.add_argument("--ld-maf005-bins", required=True)
    p.add_argument("--ld-maf005-summary", required=True)
    p.add_argument("--ld-maf001-bins", required=True)
    p.add_argument("--ld-maf001-summary", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--index-out", required=True)
    a = p.parse_args()

    reg = read_registry(a.registry)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    king = read_tsv(a.king_pairs)
    write(out / "king_kinship_vs_ibs0.html", scatter_html(
        king, "ibs0", "kinship", "IBS0 proportion", "KING-Robust kinship estimate",
        reg["king_kinship_ibs0"]["title"], reg["king_kinship_ibs0"]["caption"],
        ["sample1", "sample2", "ibs0", "kinship", "relationship"]))

    for coord, eigen, prefix in [
        (a.pca_unmasked, a.pca_unmasked_eigen, "pca_unmasked"),
        (a.pca_masked, a.pca_masked_eigen, "pca_highld"),
    ]:
        rows = read_tsv(coord)
        ev = read_eigen(eigen)
        for xpc, ypc, suffix in [(1, 2, "pc1_pc2"), (2, 3, "pc2_pc3")]:
            fid = f"{prefix}_{suffix}"
            write(out / f"{fid}.html", scatter_html(
                rows, f"PC{xpc}", f"PC{ypc}",
                f"PC{xpc} ({ev[xpc]:.2f}%)", f"PC{ypc} ({ev[ypc]:.2f}%)",
                reg[fid]["title"], reg[fid]["caption"],
                ["sample_id", "population", f"PC{xpc}", f"PC{ypc}"], "population"))

    roh = read_tsv(a.roh_individual)
    write(out / "roh_froh_ge1_5mb.html", grouped_html(
        roh, "FROH_GE1_5MB", "FROH (ROH >=1.5 Mb)",
        reg["roh_froh_ge1_5mb"]["title"], reg["roh_froh_ge1_5mb"]["caption"]))
    write(out / "roh_froh_ge5mb.html", grouped_html(
        roh, "FROH_GE5MB", "FROH from ROH >=5 Mb",
        reg["roh_froh_ge5mb"]["title"], reg["roh_froh_ge5mb"]["caption"]))
    write(out / "roh_nroh_vs_total.html", scatter_html(
        roh, "N_ROH_GE1_5MB", "TOTAL_ROH_MB_GE1_5MB",
        "Number of ROH >=1.5 Mb", "Total ROH length >=1.5 Mb (Mb)",
        reg["roh_nroh_vs_total"]["title"], reg["roh_nroh_vs_total"]["caption"],
        ["IID", "population", "N_ROH_GE1_5MB", "TOTAL_ROH_MB_GE1_5MB", "FROH_GE1_5MB", "MAX_ROH_MB"],
        "population"))

    write(out / "ld_decay_maf005.html", ld_html(
        read_tsv(a.ld_maf005_bins), read_tsv(a.ld_maf005_summary),
        reg["ld_decay_maf005"]["title"], reg["ld_decay_maf005"]["caption"]))
    write(out / "ld_decay_maf001.html", ld_html(
        read_tsv(a.ld_maf001_bins), read_tsv(a.ld_maf001_summary),
        reg["ld_decay_maf001"]["title"], reg["ld_decay_maf001"]["caption"]))

    index_rows = []
    for _, r in reg.items():
        rel_html = os.path.relpath(r["interactive_html"], Path(a.index_out).parent)
        links = [f'<a href="{esc(rel_html)}">interactive HTML</a>']
        if r.get("static_png"):
            links.append(f'<a href="{esc(os.path.relpath(r["static_png"], Path(a.index_out).parent))}">PNG</a>')
        if r.get("static_pdf"):
            links.append(f'<a href="{esc(os.path.relpath(r["static_pdf"], Path(a.index_out).parent))}">PDF</a>')
        index_rows.append(
            f"<tr><td>{esc(r['analysis'])}</td><td>{esc(r['title'])}</td>"
            f"<td>{' | '.join(links)}</td><td>{esc(r['caption'])}</td></tr>"
        )
    index = (
        '<!doctype html><html><head><meta charset="utf-8"><title>5terre figure index</title>'
        '<style>body{font-family:Arial,sans-serif;margin:24px}table{border-collapse:collapse;width:100%}'
        'th,td{border:1px solid #ccc;padding:8px;vertical-align:top}th{background:#f3f3f3}'
        'td:nth-child(4){max-width:620px}</style></head><body>'
        '<h1>5terre analysis figures</h1><p>Static manuscript-oriented figures and self-contained interactive inspection versions.</p>'
        '<table><thead><tr><th>Analysis</th><th>Figure</th><th>Files</th><th>Caption</th></tr></thead><tbody>'
        + "".join(index_rows) + '</tbody></table></body></html>'
    )
    write(a.index_out, index)

if __name__ == "__main__":
    main()
