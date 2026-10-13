# Who are our Atacama *Hyalella*? A cox1 barcode comparison

*Pipeline run 2026-10-02. Tree: `figure/07_tree_cox1.pdf` (zoomable) / `.png`.
Per-sample table: `08_nearest_relatives.tsv`; grouped: `08_nearest_relatives_grouped.tsv`.
Slides: `summary.html`.*

## The short answer

Our 40 cox1 sequences fall into three very different groups. Most of the PCR
samples (18 of 30) belong to an **unnamed high-Andean lineage** that has
previously been found in the Argentine Puna, in the Antofagasta region of
Chile and in Peru. That is the biogeographically expected answer for Atacama
amphipods. Three more PCR samples sit with the **Lake Titicaca species flock**.

The surprise is the rest: **all 10 genome-derived sequences and 9 PCR
samples are almost identical to North American lineages** (*Hyalella azteca*
and a species from a Nevada thermal spring), including lineages kept as
laboratory cultures. That is not what native Atacama animals should look
like. Before any biological interpretation, this needs checking for
contamination or sample mix-up.

None of our samples looks like a new lineage on the barcode evidence alone:
every one is at least 97 % identical to a sequence already in the public
databases.

## What we did, in plain terms

cox1 (the "DNA barcode" gene used to tell animal species apart) was compared
between our samples and every *Hyalella* cox1 we could find in the two big
public databases: NCBI GenBank (1,021 records) and BOLD, the Barcode of Life
Data System (2,088 additional records not in GenBank). Together that is
3,109 public *Hyalella* records, 2,409 of which say which country they come
from. We added one sequence each from 10 other amphipod genera to anchor
("root") the tree.

Quality control removed 8 sequences shorter than 300 bases and 58 public
records that were less than 70 % identical to our *Hyalella* reference (most
likely mis-labelled or a different gene region). All 40 of our sequences
passed. Identical sequences were collapsed, leaving 1,949 distinct sequences.

Two different comparisons were then made:

1. **Nearest relative (the table).** Each of our sequences was compared
   directly with all 1,949 public sequences, and the closest one was
   recorded with its percentage identity. This is the most direct answer to
   "what is it closest to?".
2. **Family tree (the figure).** A tree of all 3,000 sequences would be
   unreadable, so the tree uses a representative subset: one sequence per
   group of public sequences that are ≥ 95 % identical, an extra one per
   South American country, the two nearest public relatives of each of our
   samples, all our samples and the outgroups (236 sequences). The tree
   shows how the groups relate to each other and how confident we can be
   about each grouping.

### How to read the numbers

- **% identity**: the share of DNA letters that match. Within one species,
  cox1 typically differs by 0–2 %; between species, usually by more than
  3 %. We call ≥ 98 % "matches a known lineage", 96–98 % "borderline", and
  < 96 % "candidate new lineage".
- **Bootstrap support** (black dots on the tree): how often a grouping
  reappears when the data are resampled. ≥ 95 % is strong.

## Findings

### 1. Most PCR samples belong to an Andean lineage that also lives in Chile

18 PCR samples (anvi, biri, brsk, chol, deep, fern, gaze, isle, jolt, moor,
nova, oslo, quay, silt, tarn, vale, yawn, zest) form one tight group in the
tree. It includes public sequences from Argentina, Chile and Peru, and the
grouping is strongly supported (bootstrap 99 %).

- Their closest public relatives are labelled ***Hyalella* sp. O** (BOLD
  cluster BOLD:ACL3635), from **Calalaste** and **Quebrada del Diablo**,
  Argentina, at **97.2–99.4 % identity**.
- The same cluster contains unnamed *Hyalella* from **Antofagasta, Chile**
  (a lake recorded in BOLD as "Menique lagoon"). Our samples are 97.2–98.7 %
  identical to these Chilean animals. One sample (biri) is closest among
  Chilean records to an animal from the Río Salvador at Chacance in the
  Atacama Desert (96.5 %).
- Peruvian specimens labelled "3TK10", from a study of *Hyalella* speciation
  across the Andean Altiplano, also fall in this cluster.
- One sample, **tarn**, technically has its single best hit to *Hyalella*
  sp. C2 (a Titicaca lineage, 98.6 %). However, it is almost as close to
  *H.* sp. O (98.3 %) and sits inside this group in the tree.

**Meaning:** these animals are part of a widespread high-Andean / Altiplano
lineage that has been found before, but has **no formal species name**.
"*Hyalella* sp. O" is a placeholder. Our samples are not new to science
genetically, but they may belong to an undescribed species.

### 2. Three PCR samples are part of the Lake Titicaca species flock

alba, loch and wand are 98.1–99.1 % identical to *Hyalella nefrens* and
*Hyalella* sp. A1. These records come from studies of the Lake Titicaca
species flock; the GenBank records carry no country, which we took from the
publication titles. BOLD records from the Bolivian side of the lake are
similarly close (98.1–98.9 %).

Two caveats:

- The three form their own small group in the tree, and *where exactly*
  that group attaches among the Titicaca lineages is not well supported.
- The barcode does not separate Titicaca species well. The BOLD cluster
  containing *H. nefrens* also contains *H. cuprea* and *H. montforti*, and
  one of the source papers is titled "…perspectives and drawbacks of
  DNA-based identification".

So "a Titicaca-flock lineage" is solid; a species name is not.

### 3. Nineteen sequences point to North America, which needs checking

| Our samples | Nearest public lineage | Where it was found | % identity |
| --- | --- | --- | --- |
| all 10 genome-derived (ata01–ata10) | *H. azteca* complex | eastern North America; relatives with locality: Ontario, Ohio, Alberta | 99.4–99.8 |
| crsp, easy (PCR) | *H. azteca* | a laboratory culture in Lyon, France (toxicity-test study); also Germany | 98.3–99.2 |
| kelp, rove, drft (PCR) | *H. azteca* | Alberta, Canada | 98.3–99.1 |
| hope, pine, urge, xena (PCR) | *Hyalella* sp. 1 AGC-2024c | Hot Creek Spring, Nye County, Nevada, USA (a recently described thermal-spring species) | 97.4–99.0 |

In the tree these matches are well supported: 99–100 % bootstrap for the Nevada, Canadian and French groups; 81–95 % for the genome-derived sequences and their North American relatives.

**Why this is surprising:** *H. azteca* is a North American species complex.
It is also the standard amphipod kept in ecotoxicology laboratories
worldwide. The Nevada lineage is a narrow-range spring endemic. Finding
near-identical copies of all of these in Atacama samples, at the same
gene-level similarity normally seen *within* one population, is not a
plausible natural pattern without a specific explanation.

Possible explanations, roughly from most to least likely:

1. **Contamination**: DNA or animals from *H. azteca* cultures or reference
   material, in the field, the lab or the sequencing pipeline.
2. **Sample or label mix-up** between projects.
3. **Human-mediated introduction** of North American *Hyalella*. This is
   possible in principle, but would be a notable finding needing independent
   confirmation.

It is striking that **every genome-derived sequence** is North American,
while most PCR samples are Andean. If the genome assemblies were meant to
come from the same populations as the PCR samples, that alone is a red flag.
Suggested checks:

- re-extract and re-amplify a few of these samples;
- compare cox1 from the genome assemblies with Sanger PCR products from the
  same individuals;
- check lab records for *H. azteca* cultures handled alongside these samples;
- look at the morphology of the vouchers.

### 4. Are any of our samples new lineages?

On barcode evidence, no. Every sample is ≥ 97.2 % identical to a public
sequence. 36 of 40 are ≥ 98 %. The 4 borderline samples (biri 97.2 %,
hope 97.4 %, xena 97.7 %, jolt 97.9 %) are still within the range usually
seen inside one species. However:

- The main Atacama group (finding 1) is a **lineage without a species
  name**, so a formal description may still be warranted. Morphology and
  more genes would be needed.
- The Titicaca-flock samples cannot be named to species with cox1 alone.

## How confident are we?

| Statement | Confidence | Why |
| --- | --- | --- |
| 18 PCR samples belong to the high-Andean *H.* sp. O lineage found in Argentina, Chile and Peru | **High** | 97–99 % identity to many records; strongly supported group |
| 3 PCR samples belong to the Lake Titicaca flock | **High** for the flock, **low** for a species name | cox1 doesn't separate the flock's species |
| 19 sequences match North American lineages | **High** that the sequences match; **open** whether the animals really are from the Atacama | identities 97.4–99.8 %, bootstrap 81–100 %; biologically implausible |
| No clearly new species | **Moderate** | one gene; unnamed lineages exist |

General limitations:

- **One mitochondrial gene.** cox1 tells us about maternal lineages only.
- **~640 bases compared.** The genome-derived sequences are full-length
  (~1,560 bases), but the tree uses only the barcode region (638 aligned
  positions) that most public sequences cover.
- **Public metadata gaps.** About a quarter of public records give no
  country. Locality names are reported as given by the submitters (some are
  misspelled, e.g. "Menique lagoon"). Where the country was missing, we
  only *suggest* a region from the publication title, and say so.
- **The tree is a subset.** The nearest-relative table always uses all
  1,949 sequences.
- **Tree constraint.** cox1 changes so much between distant amphipod
  families that, unconstrained, some outgroups were pulled inside *Hyalella*
  (a known artefact). We therefore told the tree program that *Hyalella* is
  one group separate from the outgroups. Relationships *within* *Hyalella*
  were not constrained.

## Methods (brief)

- **Data:** NCBI E-utilities (direct HTTP; GenBank flat files for
  metadata) and the BOLD v5 portal API, downloaded 2026-10-02. BOLD records
  already in GenBank were not duplicated.
- **QC:** length ≥ 300 bp; identity to the *Hyalella* cox1 reference ≥ 70 %
  (vsearch global alignment, both strands; one public sequence was
  reverse-complemented). The reference was used only for this filter and
  is in no output. Exact duplicates were removed, keeping our own sequences
  first.
- **Alignment and tree:** MAFFT L-INS-i; alignment ends trimmed where > 50 %
  of sequences had gaps (1,570 → 638 positions). IQ-TREE 3 with ModelFinder
  (best model TIM+F+I+R5), 1,000 ultrafast bootstraps and SH-aLRT; *Hyalella*
  constrained monophyletic relative to the outgroups; rooted on the
  outgroups.
- **Nearest relatives:** vsearch exhaustive global search (identity =
  matches / aligned columns, terminal gaps ignored).
- **Reproduce:** `./src/run_pipeline.sh` (see `readme.md`).
