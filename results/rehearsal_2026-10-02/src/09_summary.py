#!/usr/bin/env python3
"""Step 09: build share/summary.html, a self-contained page for projection.

Four full-screen sections: the question, the tree, a condensed table (one
row per nearest-relative species, from 08_nearest_relatives_grouped.tsv),
five headline findings. Numbers are computed from the step 08 tables so the
page stays in sync with the pipeline. Prose lives in share/report.md.
"""
import base64
import html
import os
import sys

import pandas as pd

SHARE = os.environ["SHARE_DIR"]
PNG = os.path.join(SHARE, "figure", "07_tree_cox1.png")
FULL = os.path.join(SHARE, "08_nearest_relatives.tsv")
GROUPED = os.path.join(SHARE, "08_nearest_relatives_grouped.tsv")
META = os.path.join(os.environ["WORK_DIR"], "01_download", "01_download_metadata.tsv")
NN = os.path.join(os.environ["WORK_DIR"], "03_qc", "03_qc_nearest_public.tsv")
OUT = os.path.join(SHARE, "summary.html")

ANDEAN = {"Hyalella sp. O VI-2023a", "Hyalella nefrens", "Hyalella sp. A1",
          "Hyalella sp. C2"}
AZTECA = {"Hyalella azteca", "Hyalella sp. clade1 JL-2016"}

CSS = """
:root { --bg:#ffffff; --fg:#111111; --muted:#444; --rule:#cccccc;
        --pcr:#0062a0; --gen:#b84a00; --card:#ffffff; }
@media (prefers-color-scheme: dark) {
  :root { --bg:#0e0e0e; --fg:#f4f4f4; --muted:#c8c8c8; --rule:#555;
          --pcr:#6cc4ff; --gen:#ff9b52; --card:#ffffff; }
}
* { box-sizing:border-box; }
html { scroll-snap-type:y mandatory; }
body { margin:0; background:var(--bg); color:var(--fg);
       font-family:"Helvetica Neue",Arial,"Liberation Sans",sans-serif;
       font-size:clamp(28px,4.4vh,64px); line-height:1.25; }
section { min-height:100vh; padding:4vh 5vw; scroll-snap-align:start;
          display:flex; flex-direction:column; justify-content:center; }
h1 { font-size:clamp(44px,8vh,120px); line-height:1.1; margin:0 0 4vh; }
h2 { font-size:clamp(36px,6vh,90px); margin:0 0 2.5vh; }
p  { margin:0 0 2vh; max-width:70ch; }
.muted { color:var(--muted); }
.pcr { color:var(--pcr); font-weight:700; }
.gen { color:var(--gen); font-weight:700; }
#tree { padding:1.5vh 1vw; align-items:center; }
#tree h2 { margin-bottom:1vh; font-size:clamp(30px,4.5vh,70px); }
#tree .card { background:var(--card); border-radius:8px; padding:0.6vh;
              display:flex; justify-content:center; }
#tree img { max-height:86vh; max-width:97vw; width:auto; height:auto; display:block; }
table { border-collapse:collapse; font-size:clamp(22px,3.5vh,52px); width:100%; }
th, td { text-align:left; padding:1.1vh 1.2vw; border-bottom:2px solid var(--rule);
         vertical-align:top; }
th { font-size:0.8em; color:var(--muted); font-weight:600; }
td.num { text-align:right; white-space:nowrap; }
ul { margin:0; padding-left:1.2em; }
li { margin:0 0 2.6vh; max-width:62ch; }
a { color:var(--pcr); }
.foot { font-size:0.62em; color:var(--muted); margin-top:2vh; }
"""


def main():
    full = pd.read_csv(FULL, sep="\t", keep_default_na=False)
    grp = pd.read_csv(GROUPED, sep="\t", keep_default_na=False)
    meta = pd.read_csv(META, sep="\t", keep_default_na=False)
    img = base64.b64encode(open(PNG, "rb").read()).decode()

    n_pcr = (full.source == "pcr").sum()
    n_gen = (full.source == "genome").sum()
    hy = meta[meta.group == "hyalella"]
    n_pub = len(hy)
    pcr = full[full.source == "pcr"]
    gen = full[full.source == "genome"]
    andean = pcr[pcr.relative_species.isin(ANDEAN)]
    # Puna lineage = best hit to H. sp. O within 0.5 % of the overall top hit
    nn = pd.read_csv(NN, sep="\t", keep_default_na=False)
    nn = nn[nn.target != "*"]
    nn["organism"] = nn.target.map(meta.set_index("accession").organism)
    top = nn.groupby("query").identity.max()
    puna_best = nn[nn.organism == "Hyalella sp. O VI-2023a"].groupby("query").identity.max()
    is_puna = (puna_best.reindex(top.index) >= top - 0.5).fillna(False)
    titi = andean[~andean.our_sample.map(is_puna)]
    pcr_az = pcr[pcr.relative_species.isin(AZTECA)]
    pcr_nv = pcr[pcr.relative_species == "Hyalella sp. 1 AGC-2024c"]
    gen_az = gen[gen.relative_species.isin(AZTECA)]
    border = full[full.pct_identity < 98]

    rows = []
    for r in grp.itertuples():
        who = []
        if r.n_pcr:
            who.append(f'<span class="pcr">{r.n_pcr} PCR</span>')
        if r.n_genome:
            who.append(f'<span class="gen">{r.n_genome} genome</span>')
        where = r.closest_relative_locations.replace("*", "").replace("†", "")
        if r.n_with_chilean_hit:
            where += f"; also Chile ({r.chilean_regions})"
        pid = (f"{r.min_pct_identity:.1f}" if r.min_pct_identity == r.max_pct_identity
               else f"{r.min_pct_identity:.1f}–{r.max_pct_identity:.1f}")
        rows.append(f"<tr><td><i>{html.escape(r.relative_species)}</i></td>"
                    f"<td>{' + '.join(who)}</td><td>{html.escape(where)}</td>"
                    f'<td class="num">{pid} %</td></tr>')

    findings = [
        f"<b>{len(andean)} of {n_pcr} PCR samples are Andean <i>Hyalella</i>.</b> "
        f"Most belong to one unnamed high-Andean lineage known from the Argentine Puna "
        f"and Antofagasta, Chile; {len(titi)} sit near Lake Titicaca species.",
        f"<b>All {n_gen} genome-derived cox1 are North American <i>H. azteca</i>.</b> "
        f"They are {gen.pct_identity.min():.1f}–{gen.pct_identity.max():.1f} % identical "
        f"to haplotypes from eastern North America.",
        f"<b>{len(pcr_az) + len(pcr_nv)} PCR samples also point north.</b> "
        f"{len(pcr_az)} match <i>H. azteca</i> (Canada, a French lab culture); "
        f"{len(pcr_nv)} match a Nevada thermal-spring species.",
        f"<b>No clear new species by barcode alone.</b> Every sample is ≥ "
        f"{full.pct_identity.min():.1f} % identical to a public sequence; "
        f"{len(border)} are borderline (&lt; 98 %).",
        "<b>Before interpreting the northern matches, rule out contamination or "
        "sample mix-up.</b> Atacama animals this close to North American and "
        "lab-culture lineages is not biologically expected.",
    ]

    page = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Atacama Hyalella cox1: who are they?</title>
<style>{CSS}</style></head>
<body>
<section id="question">
  <h1>Which <i>Hyalella</i> lineages are our Atacama amphipods?</h1>
  <h2 class="muted">And where do their closest known relatives live?</h2>
  <p>We compared <span class="pcr">{n_pcr} PCR</span> and <span class="gen">{n_gen} genome</span>
  cox1 sequences with {n_pub:,} public <i>Hyalella</i> records from NCBI and BOLD.</p>
</section>
<section id="tree">
  <h2><span class="pcr">● PCR</span> &nbsp; <span class="gen">● genome</span> &nbsp;
  <span class="muted">grey = public</span></h2>
  <div class="card"><img alt="Maximum-likelihood cox1 tree of Hyalella with our samples highlighted"
  src="data:image/png;base64,{img}"></div>
</section>
<section id="table">
  <h2>Our {n_pcr + n_gen} samples, grouped by nearest known relative</h2>
  <table><thead><tr><th>Nearest public lineage</th><th>Our samples</th>
  <th>Where that relative was found</th><th>% identical</th></tr></thead>
  <tbody>{''.join(rows)}</tbody></table>
  <p class="foot">One row per sample: share/08_nearest_relatives.tsv</p>
</section>
<section id="findings">
  <h2>Five headline findings</h2>
  <ul>{''.join(f'<li>{f}</li>' for f in findings)}</ul>
  <p class="foot">Full explanation and caveats: <a href="report.md">report.md</a></p>
</section>
</body></html>
"""
    with open(OUT, "w") as fh:
        fh.write(page)
    print(f"[09] wrote {OUT} ({os.path.getsize(OUT) / 1e6:.1f} MB)")


if __name__ == "__main__":
    sys.exit(main())
