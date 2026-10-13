#!/usr/bin/env python3
"""Build share/summary.html: one self-contained page meant to be projected.

Counts come from the pipeline tables. The wording of the five findings is
interpretation: re-read it whenever the data change.
"""
import base64
import html
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cox1_common import read_tsv

groups_path, table_path, png_path, meta_path, out_path = sys.argv[1:6]
groups = read_tsv(groups_path)
samples = [r for r in read_tsv(table_path)]
meta = read_tsv(meta_path)

SOUTH_AMERICA = {'Argentina', 'Peru', 'Bolivia', 'Chile', 'Ecuador', 'Brazil', 'Colombia', 'Uruguay', 'Paraguay', 'Venezuela'}


def region(g):
    countries = set(g['nearest_public_countries'].split(', '))
    if countries & SOUTH_AMERICA:
        return 'andes'
    if 'azteca' in g['nearest_public_species'] or countries & {'USA', 'Canada'}:
        return 'north'
    return 'other'


def n_samples(gs, key=None):
    return sum(int(g['n_pcr']) + int(g['n_genome']) if key is None else int(g[key]) for g in gs)


andes = [g for g in groups if region(g) == 'andes']
north = [g for g in groups if region(g) == 'north']
n_total = len(samples)
n_placed = sum(1 for r in samples if r['percent_identity'])
min_ident = min(float(r['percent_identity']) for r in samples if r['percent_identity'])
n_chile = sum(1 for r in meta if r['set'] == 'hyalella' and r['country'] == 'Chile')
n_public = sum(1 for r in meta if r['set'] == 'hyalella')
andes_places = sorted({c for g in andes for c in g['nearest_public_countries'].split(', ')})

REGION_LABEL = {'andes': 'Andes', 'north': 'N. America', 'other': 'Unresolved'}


def row(g):
    ours = []
    if int(g['n_pcr']):
        ours.append(f'<span class="chip pcr"></span>{g["n_pcr"]} PCR')
    if int(g['n_genome']):
        ours.append(f'<span class="chip genome"></span>{g["n_genome"]} genome')
    lo, hi = g['identity_min'], g['identity_max']
    ident = f'{lo}–{hi}%' if lo != hi else f'{hi}%'
    reg = region(g)
    return (f'<tr class="{reg}"><td class="letter">{g["group"]}</td><td>{" + ".join(ours)}</td>'
            f'<td><i>{html.escape(g["nearest_public_species"].replace("Hyalella", "H."))}</i></td>'
            f'<td>{html.escape(g["nearest_public_countries"])}</td>'
            f'<td class="num">{ident}</td><td><span class="tag {reg}">{REGION_LABEL[reg]}</span></td></tr>')


findings = [
    f'<b>All {n_placed} of our {n_total} samples match lineages already in the public databases</b> (every one is at least {min_ident:.0f}% identical to a known sequence). No new lineage shows up at cox1.',
    f'<b>{n_samples(andes)} PCR samples are Andean.</b> Their closest known relatives come from the high Andes of {" and ".join(andes_places)}.',
    f'<b>{n_samples(north, "n_pcr")} PCR samples and all {n_samples(north, "n_genome")} genome samples match the North American <i>Hyalella azteca</i> group.</b> That is unexpected for the Atacama.',
    '<b>Check those North American matches before interpreting them.</b> Sample mix-up, lab-culture contamination or an introduced population are all possible.',
    f'<b>Confidence: high for which group each sample is in, low for how the groups relate.</b> Only {n_chile} of {n_public:,} public <i>Hyalella</i> cox1 records are from Chile, so "closest known" is limited by what others have sequenced.',
]

with open(png_path, 'rb') as f:
    png = base64.b64encode(f.read()).decode()

page = f'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Atacama Hyalella cox1</title>
<style>
:root {{
  color-scheme: light;
  --bg: #ffffff; --ink: #0b0b0b; --muted: #45443f; --rule: #b9b8b1; --band: #f1f0ec;
  --pcr: #c9460f; --genome: #1a57a6; --andes-bg: #fde6d8; --north-bg: #dbe8f8; --tag-ink: #0b0b0b;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    color-scheme: dark;
    --bg: #101010; --ink: #ffffff; --muted: #d2d1c8; --rule: #5c5b56; --band: #1e1e1c;
    --pcr: #ff8a55; --genome: #6fb0ff; --andes-bg: #5a2610; --north-bg: #143a6b; --tag-ink: #ffffff;
  }}
}}
:root[data-theme="dark"] {{
  color-scheme: dark;
  --bg: #101010; --ink: #ffffff; --muted: #d2d1c8; --rule: #5c5b56; --band: #1e1e1c;
  --pcr: #ff8a55; --genome: #6fb0ff; --andes-bg: #5a2610; --north-bg: #143a6b; --tag-ink: #ffffff;
}}
* {{ box-sizing: border-box; }}
html {{ font-size: clamp(20px, 2.1vw, 46px); scroll-snap-type: y proximity; }}
body {{ margin: 0; background: var(--bg); color: var(--ink);
  font-family: "Helvetica Neue", Helvetica, Arial, system-ui, sans-serif; line-height: 1.25; }}
section {{ min-height: 100vh; padding: 4vh 4vw; display: flex; flex-direction: column;
  justify-content: center; scroll-snap-align: start; border-bottom: 2px solid var(--rule); }}
.eyebrow {{ font-size: 0.8rem; font-weight: 700; letter-spacing: 0.08em; text-transform: uppercase; color: var(--muted); margin: 0 0 2vh; }}
h1 {{ font-size: 2.7rem; line-height: 1.08; margin: 0 0 4vh; text-wrap: balance; }}
h2 {{ font-size: 1.5rem; line-height: 1.1; margin: 0 0 2.5vh; }}
p {{ font-size: 1.25rem; margin: 0 0 2vh; max-width: 34em; }}
.muted {{ color: var(--muted); }}
a {{ color: inherit; text-decoration-thickness: 0.08em; text-underline-offset: 0.15em; }}
.chip {{ display: inline-block; width: 0.75em; height: 0.75em; margin-right: 0.3em; vertical-align: -0.03em; }}
.chip.pcr {{ background: var(--pcr); border-radius: 50%; }}
.chip.genome {{ background: var(--genome); }}
.key {{ font-size: 1.1rem; display: flex; gap: 2.5em; flex-wrap: wrap; font-weight: 600; }}
#tree {{ padding: 1.5vh 1vw; }}
#tree img {{ display: block; width: 100%; height: 97vh; object-fit: contain; background: #ffffff; border-radius: 6px; }}
table {{ border-collapse: collapse; width: 100%; font-size: 1.08rem; }}
th {{ text-align: left; font-size: 0.7rem; letter-spacing: 0.04em; white-space: nowrap; text-transform: uppercase; color: var(--muted);
  padding: 0 0.45em 1vh; border-bottom: 3px solid var(--ink); }}
td {{ padding: 1.25vh 0.45em; border-bottom: 1px solid var(--rule); white-space: nowrap; }}
td.letter {{ font-weight: 800; }}
td.num, th.num {{ text-align: right; font-variant-numeric: tabular-nums; }}
.tag {{ display: inline-block; padding: 0.08em 0.55em; border-radius: 0.3em; font-weight: 700; font-size: 0.9em; color: var(--tag-ink); }}
.tag.andes {{ background: var(--andes-bg); }}
.tag.north {{ background: var(--north-bg); }}
.tag.other {{ background: var(--band); }}
ol {{ margin: 0; padding-left: 1.4em; }}
li {{ font-size: 1.12rem; margin-bottom: 2.2vh; max-width: 36em; padding-left: 0.3em; }}
li::marker {{ font-weight: 800; }}
.foot {{ font-size: 0.85rem; margin-top: 2vh; }}
</style>
</head>
<body>

<section id="question">
  <div class="eyebrow">Atacama Hyalella · cox1 · preliminary v1</div>
  <h1>Which <i>Hyalella</i> lineages are our samples, and where do their closest known relatives live?</h1>
  <p>We compared one mitochondrial gene (cox1) from {n_total} of our amphipods with {n_public:,} public <i>Hyalella</i> sequences.</p>
  <div class="key">
    <span><span class="chip pcr"></span>{sum(r["source"] == "pcr" for r in samples)} PCR samples</span>
    <span><span class="chip genome"></span>{sum(r["source"] == "genome" for r in samples)} from genome assemblies</span>
  </div>
</section>

<section id="tree">
  <img alt="cox1 gene tree: our samples in orange (PCR) and blue (genome) among public Hyalella sequences, in eight lettered groups" src="data:image/png;base64,{png}">
</section>

<section id="table">
  <div class="eyebrow">Nearest public relative, by group (letters match the tree)</div>
  <table>
    <thead><tr><th>Group</th><th>Our samples</th><th>Nearest public relative</th><th>Its country</th><th class="num">DNA identity</th><th>Lineage from</th></tr></thead>
    <tbody>
    {chr(10).join("    " + row(g) for g in groups)}
    </tbody>
  </table>
  <p class="foot muted">One row per sample: <a href="table/04_nearest_public_relatives.tsv">table/04_nearest_public_relatives.tsv</a></p>
</section>

<section id="findings">
  <h2>Five findings</h2>
  <ol>
    {chr(10).join("    <li>" + f + "</li>" for f in findings)}
  </ol>
  <p class="foot muted">Full explanation and caveats: <a href="report.md">report.md</a></p>
</section>

</body>
</html>
'''
with open(out_path, 'w') as f:
    f.write(page)
print(f'[09] groups: {len(groups)} (andean {len(andes)}, north american {len(north)}); '
      f'page size {len(page) / 1e6:.1f} MB')
print(f'[09] wrote {out_path}')
