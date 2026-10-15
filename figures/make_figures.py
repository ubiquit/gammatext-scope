"""Regenerate the article figures from a reference dataset in expected/.

  python figures/make_figures.py --dataset v2 --out figures/out/v2

Produces Figure_1 (per-chromosome stop-codon difference, breast cancer, npc-DNA),
Figure_2 / Figure_3 (DNA-scope of chr1, breast cancer / melanoma, npc-DNA),
Appendix_1 (plaited DNA-scopes, npc-DNA, 3 A4 pages) and Appendix_2 (pc-DNA),
plus tables.md with the DFA summaries. Decay-dot and centromere logic mirrors
DecayGapDfa / StopCodonAlphaComparison.FindCentromereWindowRange in the C# code.
"""
import argparse, csv, os, textwrap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

DATASETS = {
    "v1": dict(bc_density="BreastCancer_Normal_vs_Tumor", bc_alpha="BreastCancer", mel_density="Melanoma_Normal_vs_Tumor", mel_alpha="Melanoma"),
    "v2": dict(bc_density="BreastCancer_v2_Normal_vs_Tumor", bc_alpha="BreastCancer_v2", mel_density="Melanoma_v3_Normal_vs_Tumor", mel_alpha="Melanoma_v3"),
}
CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]
MODE_LABEL = {"Intergenic": "npc-DNA (intergenic)", "Gene": "pc-DNA (gene bodies)", "WholeGenome": "whole DNA"}
RED, BLUE, GREY, AMBER, BAND = "#d62728", "#3c78d8", "#b8b8b8", "#f5a623", "#dbe9ff"
CLASS_COLOR = {"alfa": RED, "zipf": BLUE, "ikke": GREY}
# GRCh38 centromeres (UCSC cytoBand acen), raw bp — same table as CentromereRegions.cs
CENTROMERE = {"1": (121_700_000, 125_100_000), "2": (91_800_000, 96_000_000), "3": (87_800_000, 94_000_000), "4": (48_200_000, 51_800_000),
              "5": (46_100_000, 51_400_000), "6": (58_500_000, 62_600_000), "7": (58_100_000, 62_100_000), "8": (43_200_000, 47_200_000),
              "9": (42_200_000, 45_500_000), "10": (38_000_000, 41_600_000), "11": (51_000_000, 55_800_000), "12": (33_200_000, 37_800_000),
              "13": (0, 18_500_000), "14": (0, 18_200_000), "15": (0, 19_000_000), "16": (35_300_000, 38_400_000), "17": (22_700_000, 27_400_000),
              "18": (15_400_000, 21_500_000), "19": (24_200_000, 28_100_000), "20": (26_200_000, 30_400_000), "21": (0, 13_000_000),
              "22": (0, 17_400_000), "X": (58_100_000, 63_800_000), "Y": (10_300_000, 10_600_000)}
SLIDING, MIN_DELTA = 10, 3


# ---------------------------------------------------------------- data
def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def per_chromosome_delta(path):
    out = {}
    for r in read_tsv(path):
        out[r["Chromosome"]] = out.get(r["Chromosome"], 0) + int(r["Tumor_count"]) - int(r["Normal_count"])
    return out

def alpha_rows(exp, pair, mode, chrom):
    rows = read_tsv(os.path.join(exp, "Tumor", pair, f"StopCodonAlphaComparison_{mode}", f"StopCodonAlpha_{pair}_chr{chrom}.tsv"))
    return dict(delta=[int(r["Delta"]) for r in rows], n=[r["Normal_Class"] for r in rows], t=[r["Tumor_Class"] for r in rows],
                nstops=[int(r["Normal_StopCodons"]) for r in rows], tstops=[int(r["Tumor_StopCodons"]) for r in rows],
                start=[int(r["StartPos"]) for r in rows], length=[int(r["Length"]) for r in rows])

def decay_dots(normal, tumour):
    """Mirror of DecayGapDfa: sliding window of 10, |qualified difference| >= 3.
    Returns (dots where normal has more, dots where tumour has more)."""
    n = min(len(normal), len(tumour)); n_side, t_side = [], []
    if n < SLIDING: return n_side, t_side
    for c in range(n):
        s = max(0, c - SLIDING // 2); e = min(n, s + SLIDING)
        if e - s < SLIDING and e == n: s = max(0, e - SLIDING)
        na = sum(1 for j in range(s, e) if normal[j] == "alfa"); nz = sum(1 for j in range(s, e) if normal[j] == "zipf")
        ta = sum(1 for j in range(s, e) if tumour[j] == "alfa"); tz = sum(1 for j in range(s, e) if tumour[j] == "zipf")
        d = (na + nz) - (ta + tz)
        if abs(d) >= MIN_DELTA:
            if d > 0: n_side.append((c, "alfa" if na > nz else "zipf"))
            else: t_side.append((c, "alfa" if ta > tz else "zipf"))
    return n_side, t_side

def centromere_ranges(chrom, rows, mode="Intergenic"):
    """Centromere window ranges (inclusive index pairs), after FindCentromereWindowRange in the C# code.
    Pass 1: windows in which the tumour has no stop codons (unassembled satellite / N runs), each run expanded
    over neighbours with a tumour/normal stop-codon ratio outside 0.5–2.0. Unlike the legacy chart, which
    painted everything from the first to the last such window, each run is marked on its own, and a run that
    consists only of the chromosome's trailing partial window is ignored.
    Pass 2 (no such window): the GRCh38 cytoband position mapped proportionally — npc-DNA only."""
    ts, ns = rows["tstops"], rows["nstops"]; n = len(ts)
    zero = [i for i, v in enumerate(ts) if v == 0]
    if zero:
        runs = []
        for i in zero:
            if runs and i == runs[-1][1] + 1: runs[-1][1] = i
            else: runs.append([i, i])
        out = []
        for s, e in runs:
            if s == e == n - 1 and n > 1: continue
            while s > 0 and ns[s - 1] > 0 and not (0.5 <= ts[s - 1] / ns[s - 1] <= 2.0): s -= 1
            while e < n - 1 and ns[e + 1] > 0 and not (0.5 <= ts[e + 1] / ns[e + 1] <= 2.0): e += 1
            if out and s <= out[-1][1] + 1: out[-1] = (out[-1][0], e)
            else: out.append((s, e))
        if out: return out
    # the proportional fallback assumes intergenic coordinates (raw ≈ 1.6 × filtered); it is meaningless for
    # gene-body DNA, so pc-DNA panels mark the centromere only where the alignment itself shows it
    if mode != "Intergenic" or chrom not in CENTROMERE or not rows["start"]: return []
    total = rows["start"][-1] + rows["length"][-1]
    ratio = 1.0 / 1.6
    cs, ce = CENTROMERE[chrom][0] * ratio, CENTROMERE[chrom][1] * ratio
    margin = (ce - cs) * 0.3; cs, ce = max(0, cs - margin), min(total, ce + margin)
    idx = [i for i, (st, ln) in enumerate(zip(rows["start"], rows["length"])) if st < ce and st + ln > cs]
    return [(idx[0], idx[-1])] if idx else []


# ---------------------------------------------------------------- drawing helpers
def draw_bands(ax_t, ax_n, rows, dot_ax_t, dot_ax_n, dot_size):
    n = len(rows["delta"])
    for ax, cls in ((ax_t, rows["t"]), (ax_n, rows["n"])):
        # draw runs of equal class as single bars for speed
        i = 0
        while i < n:
            j = i
            while j + 1 < n and cls[j + 1] == cls[i]: j += 1
            ax.axvspan(i - 0.5, j + 0.5, color=CLASS_COLOR.get(cls[i], GREY), lw=0)
            i = j + 1
        ax.set_xlim(-0.5, n - 0.5); ax.set_yticks([]); ax.set_xticks([])
        for sp in ax.spines.values(): sp.set_visible(False)
    n_side, t_side = decay_dots(rows["n"], rows["t"])
    for ax, dots in ((dot_ax_t, t_side), (dot_ax_n, n_side)):
        ax.scatter([i for i, _ in dots], [0] * len(dots), s=dot_size, color=[RED if k == "alfa" else BLUE for _, k in dots], lw=0, zorder=3)
        ax.set_xlim(-0.5, n - 0.5); ax.set_ylim(-1, 1); ax.set_yticks([]); ax.set_xticks([])
        for sp in ax.spines.values(): sp.set_visible(False)

def legend_handles():
    return [Patch(color=RED, label="Alfa-qualified (Zipf + wizard-hat)"), Patch(color=BLUE, label="Zipf-qualified only"),
            Patch(color=GREY, label="unqualified"), Line2D([0], [0], color=AMBER, lw=2.5, label="centromere"),
            Line2D([0], [0], marker="o", color=BLUE, lw=0, markersize=7, label="Zipf-decay dot (3+ qualified windows ahead in 10-window region)"),
            Line2D([0], [0], marker="o", color=RED, lw=0, markersize=7, label="alfa-decay dot (same, alfa-qualified windows dominate)")]


# ---------------------------------------------------------------- Figure 1
def figure_1(exp, ds, out, pair_key="bc_density", title_sample="Breast cancer (HCC1395)", mode="Intergenic", name="Figure_1"):
    pair = ds[pair_key]
    deltas = per_chromosome_delta(os.path.join(exp, "Tumor", pair, f"StopCodonDensity_{mode}", f"{pair}_StopCodonDensity_100k_perChromosome.tsv"))
    labels = [f"chr{c}" for c in CHROMS]; values = [deltas.get(c, 0) for c in CHROMS]
    fig, ax = plt.subplots(figsize=(8.0, 9.8), dpi=200)
    y = list(range(len(labels)))[::-1]
    ax.barh(y, values, color=[RED if v > 0 else BLUE for v in values], edgecolor="#333", lw=0.5, height=0.8)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=12); ax.axvline(0, color="#444", lw=1)
    ax.set_xlabel("Stop-codon count difference (Tumour − Normal)", fontsize=13)
    ax.set_title(f"Per-chromosome stop-codon delta — {title_sample}\n{MODE_LABEL[mode]} mode", fontsize=15, weight="bold")
    ax.grid(axis="x", ls=":", alpha=0.5); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(out, f"{name}.png")); plt.close(fig)
    return dict(zip(labels, values))


def figure_1_both(exp, ds, out, mode="Intergenic", name="Figure_1_BC_MEL"):
    fig, axes = plt.subplots(1, 2, figsize=(11, 8.5), dpi=200, sharey=True)
    for ax, (key, title) in zip(axes, (("bc_density", "Breast cancer (HCC1395)"), ("mel_density", "Melanoma (COLO829)"))):
        pair = ds[key]
        deltas = per_chromosome_delta(os.path.join(exp, "Tumor", pair, f"StopCodonDensity_{mode}", f"{pair}_StopCodonDensity_100k_perChromosome.tsv"))
        values = [deltas.get(c, 0) for c in CHROMS]; y = list(range(len(CHROMS)))[::-1]
        ax.barh(y, values, color=[RED if v > 0 else BLUE for v in values], edgecolor="#333", lw=0.5, height=0.8)
        ax.set_yticks(y); ax.set_yticklabels([f"chr{c}" for c in CHROMS], fontsize=11); ax.axvline(0, color="#444", lw=1)
        ax.set_title(f"{title}" + chr(10) + f"total {sum(values):+d}", fontsize=13, weight="bold"); ax.grid(axis="x", ls=":", alpha=0.5)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        ax.set_xlabel("Stop-codon count difference (Tumour − Normal)", fontsize=11)
    fig.suptitle(f"Per-chromosome stop-codon delta — {MODE_LABEL[mode]} mode", fontsize=15, weight="bold")
    fig.tight_layout(); fig.savefig(os.path.join(out, f"{name}.png")); plt.close(fig)


# ---------------------------------------------------------------- Figure 2 / 3
def dnascope(exp, ds, out, pair_key, sample_label, pair_ids, chrom, mode, name):
    pair = ds[pair_key]; rows = alpha_rows(exp, pair, mode, chrom); n = len(rows["delta"]); d = rows["delta"]
    fig = plt.figure(figsize=(12.5, 5.9), dpi=200)
    gs = fig.add_gridspec(5, 1, height_ratios=[10, 1.1, 0.5, 1.1, 0.5], hspace=0.08, left=0.09, right=0.985, top=0.85, bottom=0.235)
    ax = fig.add_subplot(gs[0]); ax_t = fig.add_subplot(gs[1]); dot_t = fig.add_subplot(gs[2]); ax_n = fig.add_subplot(gs[3]); dot_n = fig.add_subplot(gs[4])
    lo, hi = percentile(d, 2.5), percentile(d, 97.5)
    ax.axhspan(lo, hi, color=BAND, lw=0, zorder=0)
    ylim = max(50, 4 * max(abs(lo), abs(hi)))
    clipped = [max(-ylim, min(ylim, v)) for v in d]
    ax.plot(range(n), clipped, color="black", lw=0.8, zorder=2); ax.axhline(0, color=RED, lw=1, zorder=1)
    off = [(i, v) for i, v in enumerate(d) if abs(v) > ylim]
    for i, v in off: ax.plot(i, ylim if v > 0 else -ylim, marker="^" if v > 0 else "v", color=RED, ms=5, zorder=3)
    ax.set_ylim(-ylim, ylim)
    if off:
        fig.text(0.5, 0.875, "off scale (clipped, marked ▲▼): " + ", ".join(f"window {i + 1}: {v:+d}" for i, v in off[:6]) + (" …" if len(off) > 6 else ""),
                 ha="center", fontsize=10, color=RED)
    for cen in centromere_ranges(chrom, rows, mode):
        for a in (ax, ax_t, dot_t, ax_n, dot_n): a.axvspan(cen[0] - 0.5, cen[1] + 0.5, color=AMBER, alpha=0.9, lw=0, zorder=4)
    ax.set_xlim(-0.5, n - 0.5); ax.set_xticks([]); ax.set_ylabel("Δ stop-codons\n(Tumour − Normal)", fontsize=12)
    ax.grid(axis="y", ls=":", alpha=0.5); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    draw_bands(ax_t, ax_n, rows, dot_t, dot_n, dot_size=18)
    ax_t.set_ylabel("Tumour", rotation=0, ha="right", va="center", fontsize=12, weight="bold")
    ax_n.set_ylabel("Normal", rotation=0, ha="right", va="center", fontsize=12, weight="bold")
    dot_n.set_xticks(range(0, n, 100)); dot_n.set_xticklabels([str(x) for x in range(0, n, 100)], fontsize=10); dot_n.set_xlabel("Window index", fontsize=12)
    fig.suptitle(f"Stop-codon delta and text qualification — {sample_label} chr{chrom} ({MODE_LABEL[mode]} mode)", fontsize=15, weight="bold", y=0.975)
    fig.text(0.5, 0.91, f"{pair_ids} · {n} windows · 7 000 γwords per window (≈133 000 nt)", ha="center", fontsize=11, color="#555")
    fig.legend(handles=legend_handles(), loc="lower center", ncol=3, fontsize=10, frameon=False, bbox_to_anchor=(0.5, 0.005))
    fig.savefig(os.path.join(out, f"{name}.png")); plt.close(fig)
    return n


# ---------------------------------------------------------------- Appendices
def grouped_minmax(d, size=10):
    return [(i, min(d[i:i + size]), max(d[i:i + size])) for i in range(0, len(d), size)]

def panel(fig, gs_cell, exp, ds, chrom, mode):
    """One chromosome: MEL delta + bands, BC delta + bands. Returns the off-scale note."""
    sub = gs_cell.subgridspec(10, 1, height_ratios=[6, 0.9, 0.5, 0.9, 0.5, 6, 0.9, 0.5, 0.9, 0.5], hspace=0.12)
    axes = [fig.add_subplot(sub[i]) for i in range(10)]
    notes = []
    for k, (pair_key, label) in enumerate((("mel_alpha", "MEL"), ("bc_alpha", "BC"))):
        rows = alpha_rows(exp, ds[pair_key], mode, chrom); d = rows["delta"]; n = len(d)
        a, bt, dt, bn, dn = axes[5 * k: 5 * k + 5]
        lo, hi = percentile(d, 2.5), percentile(d, 97.5)
        ylim = max(50, 4 * max(abs(lo), abs(hi)))
        a.axhspan(lo, hi, color=BAND, lw=0, zorder=0); a.axhline(0, color=RED, lw=0.7, zorder=1)
        off = []
        for i, mn, mx in grouped_minmax(d):
            cmn, cmx = max(mn, -ylim), min(mx, ylim)
            a.plot([i + 5, i + 5], [cmn, cmx], color="black", lw=1.1, solid_capstyle="butt", zorder=2)
            if mx > ylim: a.plot(i + 5, ylim, marker="^", color=RED, ms=3, zorder=3); off.append(mx)
            if mn < -ylim: a.plot(i + 5, -ylim, marker="v", color=RED, ms=3, zorder=3); off.append(mn)
        if off: notes.append(f"{label} {max(off, key=abs):+d}")
        a.set_xlim(-0.5, n - 0.5); a.set_ylim(-ylim, ylim); a.set_xticks([]); a.set_yticks([])
        a.set_ylabel(label, rotation=0, ha="right", va="center", fontsize=9, weight="bold")
        for sp in ("top", "right", "bottom"): a.spines[sp].set_visible(False)
        a.spines["left"].set_color("#aaa")
        for cen in centromere_ranges(chrom, rows, mode):
            for ax_ in (a, bt, dt, bn, dn): ax_.axvspan(cen[0] - 0.5, cen[1] + 0.5, color=AMBER, alpha=0.9, lw=0, zorder=4)
        draw_bands(bt, bn, rows, dt, dn, dot_size=6)
        bt.set_ylabel("T", rotation=0, ha="right", va="center", fontsize=7); bn.set_ylabel("N", rotation=0, ha="right", va="center", fontsize=7)
    return notes, n

def appendix(exp, ds, out, mode, name):
    pages = [CHROMS[0:8], CHROMS[8:16], CHROMS[16:24]]
    for p, chroms in enumerate(pages, start=1):
        fig = plt.figure(figsize=(8.27, 11.69), dpi=200)
        gs = fig.add_gridspec(4, 2, left=0.07, right=0.98, top=0.868, bottom=0.075, hspace=0.35, wspace=0.18)
        for k, chrom in enumerate(chroms):
            cell = gs[k % 4, k // 4]
            notes, n = panel(fig, cell, exp, ds, chrom, mode)
            bb = cell.get_position(fig)
            fig.text(bb.x0 + 0.005, bb.y1 + 0.012, f"chr{chrom}", fontsize=12, weight="bold")
            fig.text(bb.x0 + 0.075, bb.y1 + 0.012, f"{n} windows", fontsize=8, color="#666")
            if notes: fig.text(bb.x1, bb.y1 + 0.012, "off scale: " + " / ".join(notes), fontsize=8, color=RED, ha="right")
        fig.suptitle(f"Stop-codon delta and text qualification — {MODE_LABEL[mode]} mode", fontsize=14, weight="bold", x=0.03, ha="left", y=0.975)
        sub = ("Melanoma (MEL, COLO829 / COLO829BL) over breast cancer (BC, HCC1395 / HCC1395BL). "
               f"{MODE_LABEL[mode]} mode, 7 000-γword windows (≈133 000 nt). Bars: min–max of Δ stop-codons (Tumour − Normal) in 10-window groups "
               "(no peak lost); shaded band: central 95 % of Δ; red marks and notes: values beyond the panel range; orange: centromere "
               "(from the alignment, or the GRCh38 cytoband position mapped proportionally in npc-DNA mode). "
               f"Page {p} of 3.")
        fig.text(0.03, 0.952, chr(10).join(textwrap.wrap(sub, 128)), fontsize=7.5, color="#444", va="top", linespacing=1.4)
        fig.legend(handles=legend_handles(), loc="lower center", ncol=3, fontsize=7.5, frameon=False, bbox_to_anchor=(0.5, 0.008))
        fig.savefig(os.path.join(out, f"{name}_page{p}.png")); plt.close(fig)


# ---------------------------------------------------------------- tables
def dfa_table(exp, ds, out):
    lines = ["# DFA summaries (decay-marker gap series)\n", "Category rows: dots on the normal side (tumour lost qualification) / tumour side; "
             "DFA α only where the series has ≥ 30 gaps (≥ 15 with a shortened threshold); n.c. otherwise.\n"]
    for cancer, key in (("Breast cancer", "bc_alpha"), ("Melanoma", "mel_alpha")):
        pair = ds[key]
        for mode in ("Intergenic", "Gene", "WholeGenome", "WholeGenomeMerged"):
            path = os.path.join(exp, "Tumor", pair, f"DecayGapDFA_{mode}", f"DecayGapDFA_{pair}_Summary.txt")
            if not os.path.exists(path): continue
            lines.append(f"\n## {cancer} — {mode}\n\n| series | dots | gaps | mean gap | DFA α |\n|---|---|---|---|---|")
            for r in read_tsv(path):
                if not r.get("Category"): continue
                alpha = r["DFA_Alpha"]; alpha = alpha if alpha.replace(".", "", 1).isdigit() else "n.c."
                lines.append(f"| {r['Category']} | {r['DotCount']} | {r['GapCount']} | {r['MeanGap']} | {alpha} |")
    with open(os.path.join(out, "tables.md"), "w", encoding="utf-8") as fh: fh.write("\n".join(lines) + "\n")


def percentile(values, p):
    s = sorted(values); k = (len(s) - 1) * p / 100; f = int(k); c = min(f + 1, len(s) - 1)
    return s[f] + (s[c] - s[f]) * (k - f)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--dataset", default="v2"); ap.add_argument("--expected", default=os.path.join(os.path.dirname(__file__), "..", "expected"))
    ap.add_argument("--out", default=None); args = ap.parse_args()
    ds = DATASETS[args.dataset]; exp = os.path.join(args.expected, args.dataset); out = args.out or os.path.join(os.path.dirname(__file__), "out", args.dataset)
    os.makedirs(out, exist_ok=True)
    figure_1_both(exp, ds, out)
    f1 = figure_1(exp, ds, out); print("Figure_1:", {k: v for k, v in f1.items() if abs(v) >= 100})
    print("Figure_2 windows:", dnascope(exp, ds, out, "bc_alpha", "Breast cancer", "HCC1395 / HCC1395BL", "1", "Intergenic", "Figure_2"))
    print("Figure_3 windows:", dnascope(exp, ds, out, "mel_alpha", "Melanoma", "COLO829 / COLO829BL", "1", "Intergenic", "Figure_3"))
    appendix(exp, ds, out, "Intergenic", "Appendix_1"); appendix(exp, ds, out, "Gene", "Appendix_2")
    dfa_table(exp, ds, out)
    print("written to", os.path.abspath(out))

if __name__ == "__main__":
    main()
