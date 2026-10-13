# Atacama Hyalella cox1 phylogeny

Which *Hyalella* lineages are our Atacama amphipods, and where do their
closest known relatives come from? Our cox1 (30 Sanger PCR products, 10
extracted from genome assemblies) is placed among public *Hyalella* cox1 from
NCBI GenBank and BOLD, with 10 other amphipod genera as outgroups.

Conventions: `docs/conventions.md`. Paths and parameters: `config.sh`.

## setup

```bash
mamba env create -f environment.yml      # env name: hyalella_dryrun
./src/run_pipeline.sh                    # all steps; or ./src/run_pipeline.sh 6 9
```

Optional: `export NCBI_API_KEY=...` and `NCBI_EMAIL=...` before step 01.

## pipeline

| step | script | what | output |
| ---- | ------ | ---- | ------ |
| 01 | `01_download.sh/.py` | GenBank (E-utilities, GenBank flat → metadata) + BOLD v5 portal API: all *Hyalella* COI; 1 COI record each for 10 outgroup genera. BOLD records already in GenBank are not duplicated; their BOLD metadata fills GenBank gaps | `work/01_download/` |
| 02 | `02_merge.sh` | concatenate public + `data/pcr` + `data/genomes_extracted` (reference excluded) | `work/02_merge/` (sealed) |
| 03 | `03_qc.sh/.py` | drop < `MIN_LEN` bp; vsearch vs `data/reference` (both strands, reverse-complement minus-strand hits), drop < `MIN_IDENTITY` %; headers trimmed to leading ID and made unique; identical sequences deduplicated (ours kept preferentially); exhaustive nearest-public-relative search; tree set = public clustered at `CLUSTER_ID` (+1 per extra South American country) + top `TOP_HITS_PER_SAMPLE` hits per sample + all ours + outgroups | `work/03_qc/` (sealed) |
| 04 | `04_align.sh` | MAFFT L-INS-i | `work/04_align/` (sealed) |
| 05 | `05_trim.sh/.py` | trim alignment ends: drop end columns with > `MAX_GAP_FRACTION` gaps | `work/05_trim/` (sealed) |
| 06 | `06_tree.sh` | IQ-TREE 3, ModelFinder, 1000 UFBoot + 1000 SH-aLRT, constraint: *Hyalella* monophyletic vs outgroups | `work/06_tree/` |
| 08 | `08_table.sh/.py` | nearest public relative per sample (all QC-passed public seqs, not only the tree set) | `share/08_nearest_relatives*.tsv` |
| 07 | `07_plot.sh/.py` | rooted, ladderized tree; colour by source; brackets name the lineage beside each group of ours | `share/figure/07_tree_cox1.{pdf,png}` |
| 09 | `09_summary.sh/.py` | self-contained projector page | `share/summary.html` |

Step 08 runs before 07 because the plot labels groups from the table.
`share/report.md` is the written interpretation (not generated).

## sources and headers

Source is inferred from the header (see conventions): GenBank `ACC.V`, BOLD
processid `XXXX123-YY`, `_pcr_` = PCR, `genome_` = genome-extracted.
