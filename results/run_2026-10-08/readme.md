# atacama hyalella cox1 (preliminary v1)

Which *Hyalella* lineages are our Atacama samples, and where are their
closest known relatives from? cox1 from 30 PCR samples and 10 genome
assemblies, compared with public NCBI sequences.

## run

```bash
mamba env create -f environment.yml   # once; env name: hyalella-cox1-v1
./src/run_all.sh                      # ~1 min, logs to logs/
```

Parameters and paths live in `config.sh`. Conventions: `docs/conventions.md`.

## pipeline

| script | does | main output |
| ------ | ---- | ----------- |
| `src/01_download_public.sh` | NCBI E-utilities (direct HTTP): all *Hyalella* cox1 + 40 outgroup candidates, with country/locality metadata | `work/01_download/` |
| `src/02_merge.sh` | merge public + `data/pcr` + `data/genomes_extracted` (not the QC reference) | `work/02_merge/02_merge_cox1.fasta` |
| `src/03_qc.sh` | >= 300 bp; >= 70% identity to the reference (vsearch, fixes strand); dedup; unique headers | `work/03_qc/03_qc_cox1.fasta`, `03_qc_status.tsv` |
| `src/04_nearest_subset.sh` | nearest public relative of every sample (vs all public); choose tree subset | `share/table/04_nearest_public_relatives.tsv`, `work/04_subset/` |
| `src/05_align.sh` | MAFFT | `work/05_align/05_align_cox1.fasta` |
| `src/06_trim.sh` | trimAl, drop columns with > 50% gaps | `work/06_trim/06_trim_cox1.fasta` |
| `src/07_tree.sh` | FastTree GTR+Gamma | `work/07_tree/07_tree_cox1.treefile` |
| `src/08_plot.sh` | tree figure + grouped table | `share/figure/08_tree_cox1.{pdf,png}`, `share/table/08_lineage_groups.tsv` |
| `src/09_summary.sh` | projector page | `share/summary.html` |

`src/cox1_common.py` holds shared FASTA helpers and the header-pattern source
classifier (`genome_*`, `*_pcr_*`, NCBI accession).

## deliverables (`share/`)

- `figure/08_tree_cox1.pdf`, `.png` — the tree
- `table/04_nearest_public_relatives.tsv` — one row per sample of ours
- `table/08_lineage_groups.tsv` — the same, grouped by tree cluster
- `report.md` — plain-language interpretation, confidence and speed trade-offs (hand-written; update if results change)
- `summary.html` — self-contained page for projection
