#!/usr/bin/env python3
"""Step 07: plot the cox1 tree for a non-specialist audience.

- rooted on the amphipod outgroups, ladderized, rectangular
- tip labels coloured by source (public / genome-extracted / PCR); our own
  samples get a large marker, bold label and a highlight band
- big brackets on the right name the public lineage each group of our
  samples falls next to, with where its closest relatives were collected
- nodes with strong support (UFBoot >= 95) marked with a black dot
- legend uses coloured patches/markers only, at large type

Outputs: $SHARE_DIR/figure/07_tree_cox1.pdf and .png
"""
import collections
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from Bio import Phylo
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

WORK = os.environ["WORK_DIR"]
SHARE = os.environ["SHARE_DIR"]
TREE = os.path.join(WORK, "06_tree", "06_tree_cox1.treefile")
META = os.path.join(WORK, "01_download", "01_download_metadata.tsv")
TREESET = os.path.join(WORK, "03_qc", "03_qc_tree_set.tsv")
TABLE = os.path.join(SHARE, "08_nearest_relatives.tsv")
OUTDIR = os.path.join(SHARE, "figure")

# Okabe-Ito, colour-blind safe
COLOURS = {"public": "#6b6b6b", "outgroup": "#6b6b6b",
           "genome": "#D55E00", "pcr": "#0072B2"}
LEGEND_NAMES = {"public": "Public database (NCBI / BOLD)",
                "genome": "Ours: cox1 from genome assemblies",
                "pcr": "Ours: cox1 from PCR (Sanger)"}
STRONG = 95


def short_sample(n):
    if "_pcr_" in n:
        return n.split("_pcr_")[0]
    return n.replace("genome_hyalella_sp_", "").replace("_cox1_extracted", "")


def tip_label(name, src, meta, represents):
    if src in ("pcr", "genome"):
        extra = f"  (+{len(represents)} identical)" if represents else ""
        return f"{short_sample(name)}  [{'genome' if src == 'genome' else 'PCR'}]{extra}"
    m = meta.loc[name] if name in meta.index else None
    sp = m.organism if m is not None else name
    sp = "Hyalella sp. (unnamed)" if sp == "Hyalella" else sp
    where = ""
    if m is not None and m.country:
        loc = m.locality.split(",")[0].strip() if m.locality else ""
        where = f" | {m.country}" + (f", {loc}" if loc else "")
    tag = "  [outgroup]" if src == "outgroup" else ""
    return f"{sp}{where} | {name}{tag}"


def support(clade):
    """IQ-TREE labels internal nodes 'SHaLRT/UFBoot'."""
    lab = clade.name or ""
    if clade.confidence is not None and not lab:
        lab = str(clade.confidence)
    m = re.match(r"^([\d.]+)/([\d.]+)$", lab)
    if m:
        return float(m.group(2))
    return None


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    meta = pd.read_csv(META, sep="\t", keep_default_na=False).set_index("accession")
    ts = pd.read_csv(TREESET, sep="\t", keep_default_na=False).set_index("name")
    table = pd.read_csv(TABLE, sep="\t", keep_default_na=False).set_index("our_sample")

    tree = Phylo.read(TREE, "newick")
    # root on the branch separating Hyalella (+ ours) from all outgroups:
    # first root on any outgroup so the ingroup becomes a clade, then re-root on it
    outgroups = [t for t in tree.get_terminals() if ts.loc[t.name, "source"] == "outgroup"]
    ingroup = [t for t in tree.get_terminals() if ts.loc[t.name, "source"] != "outgroup"]
    tree.root_with_outgroup(outgroups[0])
    ing = tree.common_ancestor(ingroup)
    if len(ing.get_terminals()) != len(ingroup):
        print("[07] WARNING Hyalella not monophyletic; rooted on first outgroup")
    else:
        tree.root_with_outgroup(ing)
    tree.ladderize(reverse=True)

    tips = tree.get_terminals()
    n = len(tips)
    y = {t: i for i, t in enumerate(tips)}
    x = tree.depths()
    if not max(x.values()):
        x = tree.depths(unit_branch_lengths=True)

    def set_y(clade):
        if clade.is_terminal():
            return y[clade]
        ys = [set_y(c) for c in clade.clades]
        y[clade] = (min(ys) + max(ys)) / 2
        return y[clade]
    set_y(tree.root)
    xmax = max(x.values())

    # figure geometry: landscape (~16:10) so it fills a projector screen.
    # Tip labels are small (read them in the PDF); the brackets and legend
    # are the large-type layer meant to be read from the back of a room.
    H = max(12, n * 0.068 + 2.8)
    W = 30
    XR = 2.6                       # x-axis extent in units of tree depth
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0.01, 0.015, 0.98, 1 - 2.6 / H])
    ax.set_ylim(n + 4.5, -1.5)
    label_x0 = xmax * 1.02
    bx = xmax * 1.30               # bracket position, right of the tip labels
    ax.set_xlim(-xmax * 0.01, xmax * XR)
    ax.axis("off")

    # our-sample highlight bands (drawn first, under everything)
    for t in tips:
        src = ts.loc[t.name, "source"]
        if src in ("pcr", "genome"):
            ax.axhspan(y[t] - 0.5, y[t] + 0.5, xmin=0, xmax=(bx * 1.005 + xmax * 0.01) / (xmax * (XR + 0.01)),
                       color=COLOURS[src], alpha=0.13, lw=0, zorder=0)

    # branches
    for clade in tree.find_clades(order="preorder"):
        for child in clade.clades:
            ax.plot([x[clade], x[clade]], [y[clade], y[child]], color="#222", lw=0.8,
                    solid_capstyle="butt", zorder=1)
            ax.plot([x[clade], x[child]], [y[child], y[child]], color="#222", lw=0.8,
                    solid_capstyle="butt", zorder=1)
        if not clade.is_terminal():
            s = support(clade)
            if s is not None and s >= STRONG and len(clade.get_terminals()) > 1:
                ax.plot(x[clade], y[clade], "o", ms=3.2, color="black", zorder=3)

    # tips
    for t in tips:
        src = ts.loc[t.name, "source"]
        rep = [r for r in ts.loc[t.name, "represents"].split(";") if r]
        lab = tip_label(t.name, src, meta, rep)
        ours = src in ("pcr", "genome")
        ax.plot([x[t], label_x0 * 0.995], [y[t], y[t]], color="#bbbbbb", lw=0.4,
                ls=(0, (1, 2)), zorder=1)
        if ours:
            ax.plot(x[t], y[t], "o", ms=9, color=COLOURS[src], mec="white", mew=0.8,
                    zorder=4)
        ax.text(label_x0, y[t], lab, va="center", ha="left",
                fontsize=5.4 if ours else 4.4, color=COLOURS[src],
                fontweight="bold" if ours else "normal",
                fontstyle="italic" if src == "outgroup" else "normal", zorder=4)

    # brackets for runs of our samples that sit together in the tree
    # (same source, at most MAX_RUN_GAP public tips between consecutive samples)
    MAX_RUN_GAP = 2
    ours_tips = [t for t in tips if ts.loc[t.name, "source"] in ("pcr", "genome")]
    segments = []
    for t in ours_tips:
        prev = segments[-1]["tips"][-1] if segments else None
        if prev is not None and y[t] - y[prev] <= MAX_RUN_GAP + 1 \
                and ts.loc[t.name, "source"] == ts.loc[prev.name, "source"]:
            segments[-1]["tips"].append(t)
        else:
            segments.append({"tips": [t]})
    for seg in segments:
        species = table.loc[[t.name for t in seg["tips"]], "relative_species"]
        names = [s.replace("Hyalella ", "H. ") for s in species.value_counts().index]
        seg["group"] = names[0] + (" & others" if len(names) > 1 else "")
    label_slots = []
    for seg in segments:
        y0, y1 = y[seg["tips"][0]], y[seg["tips"][-1]]
        rows = table.loc[[t.name for t in seg["tips"]]]
        locs = []
        for r in rows.itertuples():
            if r.relative_country != "not recorded":
                c = r.relative_country
                loc = r.relative_locality.split(",")[0] if r.relative_locality != "not recorded" else ""
                locs.append(f"{c}, {loc}" if loc else c)
            elif r.relative_context_from_study:
                locs.append(r.relative_context_from_study.replace(" (from publication title)", ""))
            else:
                locs.append("location not recorded")
            # a Chilean record almost as close as the top hit is worth naming
            if r.closest_chilean_pct_identity != "" and \
                    float(r.closest_chilean_pct_identity) >= r.pct_identity - 1.0:
                locs.append("Chile, " + r.closest_chilean_location.split(",")[0])
        locs = [loc for loc, _ in collections.Counter(locs).most_common()]
        where = "; ".join(locs[:2]) + (" …" if len(locs) > 2 else "")
        pid = f"{rows.pct_identity.min():.1f}–{rows.pct_identity.max():.1f} %" \
            if rows.pct_identity.min() != rows.pct_identity.max() \
            else f"{rows.pct_identity.min():.1f} %"
        label_slots.append({"y0": y0, "y1": y1, "yc": (y0 + y1) / 2,
                            "title": f"{len(seg['tips'])} × near {seg['group']}",
                            "sub": f"{where}\n{pid} identical",
                            "colour": COLOURS[ts.loc[seg['tips'][0].name, 'source']]})
    # repel labels vertically so large text never overlaps
    min_gap = 27.0  # in tip units (~1.8 in of label block)
    for i, s in enumerate(label_slots):
        s["ly"] = s["yc"] if i == 0 else max(s["yc"], label_slots[i - 1]["ly"] + min_gap)
    if label_slots and label_slots[-1]["ly"] > n - 3:   # pushed off the bottom
        label_slots[-1]["ly"] = n - 3
        for i in range(len(label_slots) - 2, -1, -1):
            label_slots[i]["ly"] = min(label_slots[i]["ly"],
                                       label_slots[i + 1]["ly"] - min_gap)
    for s in label_slots:
        c = s["colour"]
        ax.plot([bx, bx], [s["y0"] - 0.4, s["y1"] + 0.4], color=c, lw=5,
                solid_capstyle="butt")
        ax.plot([bx, bx + xmax * 0.06], [s["yc"], s["ly"]], color=c, lw=2.5)
        ax.text(bx + xmax * 0.07, s["ly"] - 0.9, s["title"], fontsize=40,
                fontweight="bold", color=c, va="bottom", ha="left")
        ax.text(bx + xmax * 0.07, s["ly"] + 0.9, s["sub"], fontsize=28,
                color="#222", va="top", ha="left", linespacing=1.25)

    # scale bar
    sb = 0.05
    ax.plot([0, sb], [n + 3.2, n + 3.2], color="black", lw=2.5)
    ax.text(sb * 1.3, n + 3.2, f"{sb} substitutions per site (≈ {sb:.0%} DNA difference)",
            ha="left", va="center", fontsize=20)

    # legend + reading guide (header band)
    handles = [Line2D([], [], marker="o", ls="", ms=30, color=COLOURS[k],
                      mec="white", label=LEGEND_NAMES[k]) for k in ("pcr", "genome")]
    handles.append(Patch(facecolor=COLOURS["public"], edgecolor="none",
                         label=LEGEND_NAMES["public"]))
    handles.append(Line2D([], [], marker="o", ls="", ms=13, color="black",
                          label=f"Branch point supported in ≥ {STRONG}% of bootstrap replicates"))
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.01, 0.998),
               ncol=2, fontsize=30, frameon=False, handlelength=1.4,
               handleheight=1.2, columnspacing=2.5, labelspacing=0.5)
    fig.text(0.99, 0.997,
             "How to read: each line ending on the right is one DNA sequence.\n"
             "Sequences joined by short paths are close relatives.\n"
             "Rooted on 10 other amphipods.",
             ha="right", va="top", fontsize=22, color="#222", linespacing=1.4)

    for ext, dpi in (("pdf", None), ("png", 110)):
        out = os.path.join(OUTDIR, f"07_tree_cox1.{ext}")
        fig.savefig(out, dpi=dpi, facecolor="white")
        print(f"[07] wrote {out}")
    print(f"[07] {n} tips; {len(ours_tips)} ours; {len(label_slots)} bracket groups")


if __name__ == "__main__":
    sys.exit(main())
