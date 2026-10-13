# Which *Hyalella* lineages are our Atacama samples? (cox1, preliminary v1)

## In one paragraph

We compared one mitochondrial gene (cox1, the standard animal "DNA barcode") from our 40 samples with every comparable *Hyalella* sequence in NCBI GenBank. Every sample matched a lineage that is already in the public record at 97% identity or better, so cox1 gives **no evidence of a new lineage** in this batch. The samples split into two very different sets. **21 of the 30 PCR samples are Andean**: their closest known relatives were collected in north-west Argentina and Peru. The other **9 PCR samples and all 10 genome-derived sequences match the North American *Hyalella azteca* species complex**, which is not what we would expect from the Atacama and needs to be checked before anyone interprets it.

## What we did, in plain terms

cox1 mutates fast enough that different species usually differ by several percent, while members of one species usually differ by less than 2–3%. So "how similar is our sequence to the closest one in the database, and where was that one collected?" is a quick first answer to "what is it, and where are its relatives?".

- **Our data:** 30 PCR (Sanger) sequences and 10 cox1 sequences extracted from our genome assemblies. All 40 passed quality control.
- **Public data:** 1,021 *Hyalella* cox1 records from GenBank plus 40 other amphipods as candidate outgroups (distant relatives used to orient the tree). After quality control and removal of exact duplicates, 869 unique public sequences remained.
- **Nearest relative:** each of our sequences was compared against all 869 public sequences. The table reports the best match.
- **Tree:** built from a subset of 92 sequences (our 40, their nearest public relatives, one representative of each formally named *Hyalella* species, and 10 outgroups).

## What the tree shows

Our samples fall into eight tight clusters, lettered A–H in the figure (`figure/08_tree_cox1.pdf`). Within each cluster the sequences are nearly identical to one another and to a public sequence.

| Group | Our samples | Nearest public relative | Where that relative is from | Identity |
|---|---|---|---|---|
| A | 18 PCR | *Hyalella* sp. O (an unnamed lineage, GenBank label "sp. O VI-2023a") | Argentina; two of the three records give the locality Calalaste | 97.2–99.4% |
| B | 3 PCR | *Hyalella nefrens* | Peru (no locality given) | 98.1–99.1% |
| C | 3 PCR | *Hyalella azteca* | Canada (Alberta, Manitoba) | 98.3–99.1% |
| D | 2 PCR | *Hyalella azteca* | France (Lyon) | 98.3–99.2% |
| F | 4 PCR | *Hyalella* sp. 1 / *H. azteca* | USA (Hot Creek Spring, Nevada) | 97.4–99.0% |
| E, G, H | 10 genome | *Hyalella azteca* | not recorded in GenBank | 99.4–99.8% |

The per-sample version is `table/04_nearest_public_relatives.tsv`; the grouped version with sample names is `table/08_lineage_groups.tsv`.

## Known species or new lineages?

**None of the 40 looks like a new lineage on this evidence.** The lowest identity to a public sequence is 97.2%, and most are above 98%. We used 97% as the rule-of-thumb line for "same lineage as something already known" and 90% as the line below which we would call a sample a candidate new lineage; nothing came close to the second line.

Two qualifications matter:

- **Group A matches a lineage that has a sequence but no name.** "*Hyalella* sp. O" is a placeholder used by the authors who deposited it. So 18 of our samples belong to a lineage that is *known to exist* but is, as far as GenBank shows, *not formally described*. That is still interesting; it just is not new to science at the cox1 level.
- **A 97% line is a convention, not a species test.** Samples at 97–98% (for example `biri`, `hope`, `xena`) sit at the edge. Whether they are the same species as their nearest relative needs more genes and morphology.

## Where the closest relatives come from

**The Andean set (groups A and B, 21 PCR samples).** Group A's relatives are from Argentina, with Calalaste given as the locality for two records; that places them in the high-Andean plateau region that adjoins the Atacama. Group B's relative, *H. nefrens*, is recorded from Peru. In the tree both groups sit among other species sampled in Peru, Bolivia and Ecuador. Geographically this is the sensible result: our Atacama animals belong to lineages of the central Andes.

**The North American set (groups C–H, 9 PCR + 10 genome samples).** These match *H. azteca* and its close relatives, a species complex native to North America. The matching public records come from Canada and the USA. The French record (group D) and its German and Swiss neighbours in the tree are *H. azteca* too, a species widely kept in laboratories for toxicity testing, so those are most likely laboratory stocks rather than wild European animals. The genome-derived sequences match GenBank records that carry no country at all.

**This second set should be treated as a red flag, not a finding.** All ten genome assemblies and nine PCR samples matching a North American species at 97–99.8% is hard to explain by natural distribution in the Atacama. More likely explanations, in the order we would check them:

1. sample or label mix-up with *H. azteca* material (a common lab organism);
2. contamination in the lab or in the sequencing run;
3. the cox1 "extracted" from assemblies was pulled out with, or assembled against, an *H. azteca* reference;
4. a genuinely introduced *H. azteca* population. This would be notable, and needs field confirmation before being claimed.

Note that not one genome sample falls in an Andean group, and not one of the 21 Andean PCR samples has a genome counterpart. Cross-checking which individuals went into the genome assemblies is the first thing to do.

## How confident are we?

| Statement | Confidence | Why |
|---|---|---|
| Which group (A–H) each sample belongs to | High | Identities of 97–99.8% over roughly 500–1,500 aligned bases; the tree and the direct comparison agree; cluster branches are well supported. |
| "No new lineage at cox1" | Moderate–high | Clear for this gene. One mitochondrial gene cannot rule out differences elsewhere in the genome. |
| Group A and B are central-Andean | Moderate–high | Clear nearest relatives with country data, but only a handful of public records, and few with precise localities. |
| The North American matches reflect animals living in the Atacama | Low | Biologically unexpected; several mundane explanations have not been excluded. |
| How the groups relate to each other (deep branches of the tree) | Low | cox1 saturates at this depth, and the tree was built fast. In this tree *Hyalella* does not even separate cleanly from the outgroups, which is a known limitation of the gene, not a discovery. |

**The biggest limit is the reference database, not the analysis.** Only 2 of the 1,021 public *Hyalella* cox1 records are from Chile, and about half of all records (525) carry no country. "Closest known relative from Argentina" means "closest among those someone has sequenced and deposited", not "closest that exists".

## What was traded away for speed

The whole pipeline runs in about one minute (46 s of that is the download). To get there:

- **Fast tree method.** FastTree (GTR+Γ) with its built-in quick branch-support values. No model selection, no bootstrap, no second method to cross-check. Support dots on the figure are indicative only.
- **Subset tree.** 92 of 909 quality-passed sequences are in the tree. The nearest-relative table does **not** use this shortcut; it searches all 869 public sequences.
- **Backbone of named species only.** Unnamed public lineages ("*Hyalella* sp. X") appear in the tree only when they are a nearest relative of one of our samples.
- **NCBI only.** BOLD was not queried. BOLD may hold South American records that GenBank lacks, which matters given how thin Chile is.
- **No whole mitochondrial genomes.** Public records longer than 2,000 bases were excluded rather than having cox1 cut out of them.
- **One gene.** cox1 only; no nuclear markers.
- **Simple alignment and trimming.** MAFFT automatic mode; columns with gaps in more than half the sequences removed (applied along the whole alignment, which in practice trims it to the 636-base barcode region). The alignment was not checked by eye or for reading-frame problems.

Other details worth knowing:

- Two public records failed the 70% identity check and were dropped (one *Hyalella longipes* record, one outgroup). 190 public sequences were exact duplicates of another and were collapsed. None of our 40 samples was dropped or collapsed.
- Two sequences were on the opposite strand and were reverse-complemented during QC.
- Where several public sequences tied as the nearest relative, the table reports the one with the best location data and lists the tie count.
- The QC reference sequence was used only for the identity check; it is not in the merged data, the alignment or the tree.

## Suggested next steps

1. Trace the 19 "North American" samples back through the lab records before any further analysis.
2. Query BOLD and add mitogenome-derived cox1 to strengthen the South American reference set.
3. For groups A and B, rerun with all public Andean sequences, a slower tree method with bootstrap support, and at least one nuclear marker.
