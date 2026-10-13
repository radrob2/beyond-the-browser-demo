#!/usr/bin/env python3
"""Draw the cox1 tree for a lay audience and write the grouped lineage table.

Every sample of ours is a tip, coloured by source. Public-only clades are
collapsed to one triangle each. Tight clades containing our samples get a
lettered bracket naming the nearest public relative and its country.
"""
import collections
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from Bio import Phylo

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cox1_common import classify, read_tsv

tree_path, meta_path, table_path, out_prefix, groups_path = sys.argv[1:6]
max_depth, support_shown = float(sys.argv[6]), float(sys.argv[7])

INK, MUTED, GRID = '#0b0b0b', '#52514e', '#c9c8c2'
COLOUR = {'pcr': '#d9531e', 'genome': '#1f63b8', 'public': MUTED}
MARKER = {'pcr': 'o', 'genome': 's'}
SOURCE_LABEL = {'pcr': 'Our samples: PCR (Sanger)', 'genome': 'Our samples: from genome assemblies',
                'public': 'Public database sequences (NCBI)'}

meta = {r['accession']: r for r in read_tsv(meta_path)}
table = {r['our_sample']: r for r in read_tsv(table_path)}
nearest = {r['nearest_public_relative'] for r in table.values()}

# ---- root on the largest outgroup-only clade, then ladderize ----------------
tree = Phylo.read(tree_path, 'newick')
is_out = lambda n: n in meta and meta[n]['set'] == 'outgroup'
is_own = lambda n: classify(n) != 'public'
tree.root_with_outgroup(next(t for t in tree.get_terminals() if is_own(t.name)))
out_clades = [c for c in tree.find_clades() if all(is_out(t.name) for t in c.get_terminals())]
root_clade = max(out_clades, key=lambda c: c.count_terminals())
n_root_out = root_clade.count_terminals()
tree.root_with_outgroup(root_clade)
tree.ladderize()

# ---- collapse clades with none of our samples and no nearest relative -------
def collapsible(clade):
    return not any(is_own(t.name) or t.name in nearest for t in clade.get_terminals())

depth = tree.depths()
rows = []          # leaves of the drawn tree: tips or collapsed clades


def species_summary(tips, limit=3):
    names = collections.Counter(meta[t.name]['organism'] for t in tips)
    countries = sorted({meta[t.name]['country'] for t in tips if meta[t.name]['country']})
    if all(is_out(t.name) for t in tips):
        what = f'{len(tips)} other amphipod genera (outgroup)'
    else:
        short = [n.replace('Hyalella ', 'H. ') for n, _ in names.most_common(limit)]
        more = f' +{len(names) - limit} more' if len(names) > limit else ''
        what = f'{len(tips)} public: ' + ', '.join(short) + more
    where = ', '.join(countries[:4]) + ('…' if len(countries) > 4 else '')
    return what + (f' · {where}' if where else '')


def layout(clade):
    """Assign y to every drawn node; return y."""
    if clade.is_terminal() or (collapsible(clade) and clade is not tree.root):
        clade.y = len(rows)
        clade.drawn_leaf = True
        rows.append(clade)
        return clade.y
    clade.drawn_leaf = False
    ys = [layout(c) for c in clade.clades]
    clade.y = (min(ys) + max(ys)) / 2
    return clade.y


layout(tree.root)
n_rows = len(rows)

# ---- lettered groups: tight clades that contain our samples -----------------
def clade_height(c):
    return max(depth[t] for t in c.get_terminals()) - depth[c]


groups = []


def find_groups(clade):
    tips = clade.get_terminals()
    if not any(is_own(t.name) for t in tips):
        return
    if clade.is_terminal() or clade_height(clade) <= max_depth:
        groups.append(clade)
        return
    for c in clade.clades:
        find_groups(c)


find_groups(tree.root)
groups.sort(key=lambda c: min(t.y for t in c.get_terminals() if hasattr(t, 'y')))

group_rows = []
for i, g in enumerate(groups):
    own = [t.name for t in g.get_terminals() if is_own(t.name)]
    recs = [table[s] for s in own]
    species = collections.Counter(r['relative_species'] for r in recs).most_common(1)[0][0]
    countries = collections.Counter(r['relative_country'] for r in recs if r['relative_country'] != 'not recorded')
    idents = [float(r['percent_identity']) for r in recs]
    g.letter = chr(ord('A') + i)
    g.n_pcr = sum(classify(s) == 'pcr' for s in own)
    g.n_genome = sum(classify(s) == 'genome' for s in own)
    g.species = species
    g.where = ', '.join(c for c, _ in countries.most_common()) or 'country unknown'
    g.ident = f'{min(idents):.1f}–{max(idents):.1f}%' if min(idents) != max(idents) else f'{max(idents):.1f}%'
    group_rows.append({
        'group': g.letter, 'n_pcr': g.n_pcr, 'n_genome': g.n_genome,
        'nearest_public_species': species, 'nearest_public_countries': g.where,
        'identity_min': f'{min(idents):.1f}', 'identity_max': f'{max(idents):.1f}',
        'our_samples': ', '.join(sorted(own))})
with open(groups_path, 'w') as f:
    cols = list(group_rows[0])
    f.write('\t'.join(cols) + '\n')
    for r in group_rows:
        f.write('\t'.join(str(r[c]) for c in cols) + '\n')

# ---- draw -------------------------------------------------------------------
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Helvetica Neue', 'Helvetica', 'Arial', 'DejaVu Sans'],
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})
W, H = 20, 11.25                                  # 16:9, fills a projector
fig = plt.figure(figsize=(W, H), facecolor='white')
ax = fig.add_axes([0.012, 0.035, 0.976, 0.80])
xmax = max(depth[c] for c in rows)
ax.set_xlim(-0.01 * xmax, xmax * 2.50)            # right 60 % holds labels + brackets
ax.set_ylim(n_rows - 0.3, -0.7)
ax.axis('off')
tip_fs = min(9.5, 0.80 * 72 * (H * 0.80) / n_rows)
LW = 1.3


def draw(clade):
    x = depth[clade]
    if clade.drawn_leaf:
        return
    ys = [c.y for c in clade.clades]
    ax.plot([x, x], [min(ys), max(ys)], color=INK, lw=LW, solid_capstyle='round', zorder=2)
    for c in clade.clades:
        ax.plot([x, depth[c]], [c.y, c.y], color=INK, lw=LW, solid_capstyle='butt', zorder=2)
        draw(c)
    if clade is not tree.root and clade.confidence is not None and clade.confidence >= support_shown:
        ax.plot(x, clade.y, 'o', ms=4.2, color=INK, zorder=3)


draw(tree.root)

label_x = {}
for leaf in rows:
    x, y = depth[leaf], leaf.y
    if not leaf.is_terminal():                    # collapsed public clade
        far = x + clade_height(leaf)
        ax.add_patch(Polygon([[x, y], [far, y - 0.36], [far, y + 0.36]], closed=True,
                             facecolor=GRID, edgecolor=MUTED, lw=0.8, zorder=2))
        ax.text(far + 0.008 * xmax, y, species_summary(leaf.get_terminals()), va='center', ha='left',
                fontsize=tip_fs, color=MUTED, style='italic')
        continue
    name = leaf.name
    src = classify(name)
    if src == 'public':
        m = meta[name]
        place = m['country'] + (f": {m['locality']}" if m['locality'] else '') if m['country'] else 'country not recorded'
        if len(place) > 28:
            place = place[:27] + '…'
        prefix = 'Outgroup: ' if is_out(name) else ''
        text = f"{prefix}{m['organism']} · {place} ({name})"
        ax.text(x + 0.008 * xmax, y, text, va='center', ha='left', fontsize=tip_fs, color=COLOUR['public'])
    else:
        ax.plot(x, y, MARKER[src], ms=tip_fs * 0.95, color=COLOUR[src], mec='white', mew=0.6, zorder=4)
        short = name.replace('genome_hyalella_sp_', 'genome ').replace('_cox1_extracted', '')
        ax.text(x + 0.016 * xmax, y, short, va='center', ha='left', fontsize=tip_fs + 1,
                color=COLOUR[src], fontweight='bold')

# brackets with large group labels
fig.canvas.draw()
bx = xmax * 1.37
group_fs = 18
min_gap = group_fs * 1.25 / (72 * H * 0.80) * (n_rows + 0.4)     # one text line, in rows
placed = []
for g in groups:
    ys = [t.y for t in g.get_terminals() if hasattr(t, 'y')] or [g.y]
    y0, y1 = min(ys) - 0.38, max(ys) + 0.38
    colour = COLOUR['pcr'] if g.n_pcr >= g.n_genome else COLOUR['genome']
    ax.plot([bx, bx], [y0, y1], color=colour, lw=5, solid_capstyle='butt', zorder=3)
    ax.axhspan(y0, y1, xmin=0, xmax=(bx - ax.get_xlim()[0]) / (ax.get_xlim()[1] - ax.get_xlim()[0]),
               color=colour, alpha=0.07, lw=0, zorder=0)
    ty = (y0 + y1) / 2
    if placed and ty - placed[-1] < min_gap:
        ty = placed[-1] + min_gap
    placed.append(ty)
    ours = ' + '.join(p for p in (f'{g.n_pcr} PCR' if g.n_pcr else '', f'{g.n_genome} genome' if g.n_genome else '') if p)
    ax.text(bx + 0.02 * xmax, ty, f'{g.letter}', va='center', ha='left', fontsize=group_fs + 2,
            fontweight='bold', color=colour)
    ax.text(bx + 0.065 * xmax, ty, f'{ours}  ·  nearest: {g.species.replace("Hyalella", "H.")}, {g.where}  ·  {g.ident}',
            va='center', ha='left', fontsize=group_fs, color=INK)

# scale bar
sb = 0.05
ax.plot([0, sb], [n_rows - 0.9, n_rows - 0.9], color=INK, lw=2.5)
ax.text(0, n_rows - 1.5, '5% DNA difference', ha='left', va='bottom', fontsize=11, color=INK)

# title, reading guide, legend (real marker handles; no text glyph keys)
fig.text(0.012, 0.955, 'Where do our Atacama Hyalella fit?  cox1 gene tree', fontsize=27, fontweight='bold', color=INK, va='center')
fig.text(0.012, 0.912, 'How to read it: every line end is one DNA sequence. Sequences joined by short branches are close relatives; '
         'longer branches mean more DNA difference.', fontsize=14.5, color=MUTED, va='center')
handles = [
    Line2D([], [], marker='o', ls='', ms=15, color=COLOUR['pcr'], label=SOURCE_LABEL['pcr']),
    Line2D([], [], marker='s', ls='', ms=14, color=COLOUR['genome'], label=SOURCE_LABEL['genome']),
    Line2D([], [], marker='_', ls='', ms=18, mew=3, color=COLOUR['public'], label=SOURCE_LABEL['public']),
    Line2D([], [], marker='<', ls='', ms=14, color=GRID, mec=MUTED, label='Public sequences, grouped'),
    Line2D([], [], marker='o', ls='', ms=7, color=INK, label=f'Well-supported split (≥ {support_shown})'),
]
leg = fig.legend(handles=handles, loc='center left', bbox_to_anchor=(0.006, 0.868), ncol=5, frameon=False,
                 fontsize=14.5, handletextpad=0.3, columnspacing=1.6, borderaxespad=0)
for text, key in zip(leg.get_texts(), ['pcr', 'genome', 'public', None, None]):
    text.set_color(COLOUR[key] if key in ('pcr', 'genome') else INK)
    if key in ('pcr', 'genome'):
        text.set_fontweight('bold')
fig.text(0.988, 0.012, f'Right-hand labels: lettered group · how many of our samples · nearest public relative, its country · % DNA identity to it.   '
         f'FastTree GTR+Γ on {tree.count_terminals()} sequences; rooted on {n_root_out} outgroup amphipods; deep branches are unreliable.',
         fontsize=10.5, color=MUTED, ha='right', va='bottom')

fig.savefig(out_prefix + '.pdf')
fig.savefig(out_prefix + '.png', dpi=200)
print(f'[08] drawn rows: {n_rows} (from {tree.count_terminals()} tips); groups: {len(groups)}')
for r in group_rows:
    print(f"[08]   {r['group']}: {r['n_pcr']} pcr + {r['n_genome']} genome -> {r['nearest_public_species']} "
          f"[{r['nearest_public_countries']}] {r['identity_min']}-{r['identity_max']}%")
print(f'[08] wrote {out_prefix}.pdf / .png and {groups_path}')
