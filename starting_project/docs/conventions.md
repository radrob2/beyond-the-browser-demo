# project conventions

Personal conventions for bioinformatics projects. Copy relevant sections to project readmes.

Last updated: 2026-05-06

---

## directory structure

| directory | purpose | git | notes |
| --------- | ------- | --- | ----- |
| `src/` | scripts (bash, python, R) | tracked | numbered by pipeline step |
| `data/` | input data from collaborators/sequencing | ignored | raw, unmodified; **drop `.claude-sealed` marker** |
| `work/` | pipeline intermediate outputs | ignored | organized by step; mark sealed if it holds DNA |
| `work/metadata-redacted/` | redacted versions of sealed metadata | ignored | Claude-readable; rerun `redact_metadata.py` on update |
| `share/` | final deliverables (figures, tables) | tracked | publication-ready |
| `logs/` | log files | ignored | timestamped |
| `bin/` | third-party tools only | ignored | not custom scripts |
| `docs/` | documentation, presentations | tracked | project docs |

## file naming

- **lowercase only**
- pattern: `{step}_{description}_{details}.{ext}`
- examples:
  - `01_download_baetidae_coi.fasta`
  - `02_align_coi.fasta`
  - `03_tree_coi.treefile`

## fasta files

- **single-line sequences**: each record is a `>header` line followed by
  one long sequence line. No 60- or 80-char line wrapping.
- **uppercase letters only**: `ACGT` plus IUPAC ambiguity codes
  (`N`, `Y`, `R`, `S`, `W`, `K`, `M`, etc.). Lowercase from NCBI or
  collaborators gets uppercased.
- **why**: single-line + uppercase makes FASTA files greppable, sed-able,
  and consistent across NCBI / collaborator / mitogenome inputs.
- **utility**: a self-contained `normalize_fasta.py` ships at
  `utils/normalize_fasta.py` in each bioinformatics project. Idempotent.
  Usage:

  ```bash
  ./utils/normalize_fasta.py INPUT.fa [OUTPUT.fa]
  ```

  If `OUTPUT.fa` is omitted, rewrites `INPUT.fa` in place.

## headers and source attribution

- **Preserve original headers.** Don't rewrite, reformat, or re-prefix
  headers from NCBI downloads or input FASTAs in `data/`. Headers are the
  primary key downstream; pipelines should pass them through unchanged
  except for whitespace/description trimming where strictly required (e.g.,
  many tree-builders refuse spaces — keep the leading accession or
  identifier token only).
- **Source attribution comes from header pattern**, not from re-tagging.
  For amphipod-pipeline-style data:
  - NCBI accessions look like `>{accession}.{version} ...` (alphanumeric
    with a dot-version, e.g. `>KU870690.1`).
  - PCR products from `data/pcr/` have headers containing `_pcr_`
    (e.g., `>alba_pcr_2025-03-01`).
  - Genome-extracted sequences from `data/genomes_extracted/` have headers
    starting with `genome_` (e.g., `>genome_hyalella_sp_ata01_cox1_extracted`).
- **Why:** preserving headers means downstream scripts can classify
  sequences from the header alone — no sidecar metadata files, no
  reheadering step that can drift between runs.

## external data fetches

NCBI E-utilities (eutils.ncbi.nlm.nih.gov) is the primary external data
source. **Use direct HTTP via `curl` or Python's `urllib`/`requests` —
never `entrez-direct` (`esearch | efetch | xtract` shell pipeline) for
bulk queries.** The direct-HTTP route is orders of magnitude faster for
the metadata-heavy fetches typical in phylogenetics work.

Pattern for bulk fetches:

1. `esearch.fcgi?db=…&term=…&usehistory=y` — get `WebEnv` + `QueryKey`.
2. `efetch.fcgi?db=…&WebEnv=…&query_key=…&rettype=fasta&retmode=text` —
   pull all matching records via the history server in one call.

If `NCBI_API_KEY` is in the environment, append `&api_key=$NCBI_API_KEY`
to every E-utilities URL to raise the per-second rate limit.

## scope: stay within the project

When writing scripts, drawing on conventions, or hunting for reference
implementations, stay within the current project directory. Do not read,
copy from, or reference code in sibling projects under `~/projects/` —
that pattern is convenient in interactive development but harms
reproducibility (the demo audience or a collaborator who clones this
project doesn't have those sibling projects). Anything that should be
shared belongs in this conventions doc or in the project's own `utils/`.

## script conventions

1. **Numbering**: Scripts numbered by pipeline order (`01_`, `02_`, etc.)
2. **Config**: All scripts must `source config.sh` for paths and parameters
3. **Shebang**: Use `#!/bin/bash` for shell, `#!/usr/bin/env python3` for python
4. **Error handling**: Use `set -euo pipefail` in bash scripts
5. **Output**: Scripts should echo progress and final counts
6. **Paths**: Use `$WORK_DIR`, `$DATA_DIR`, `$SHARE_DIR` from config.sh

## config.sh

Central configuration file containing:

- `PROJECT_DIR`, `DATA_DIR`, `WORK_DIR`, `SHARE_DIR`, `LOGS_DIR`
- `THREADS` - number of CPU threads for parallel tools
- `CONDA_ENV` - environment name
- `activate_conda()` - function to activate conda

Source in scripts:
```bash
source "$(dirname "$0")/../config.sh"
```

## adding new scripts

1. Number sequentially based on pipeline position
2. Add `source "$(dirname "$0")/../config.sh"` at top
3. Use `set -euo pipefail` for bash
4. Create output directory with `mkdir -p "$OUTDIR"`
5. Echo progress messages
6. Update readme if adding new pipeline steps

## adding new data

1. Place raw input files in `data/` (gitignored)
2. Downloaded/processed files go in `work/{step}/`
3. Never modify files in `data/` - always copy to `work/`

## git tracking

**Tracked:**
- `src/` - all scripts
- `share/` - final figures, tables, deliverables
- `config.sh`, `environment.yml`, `readme.md`

**Ignored:**
- `data/` - large input files
- `work/` - reproducible intermediate files
- `logs/` - regeneratable logs
- `bin/` - third-party binaries

## logging

- Log files go in `logs/` directory (gitignored)
- Naming: `{script}_{YYYYMMDD}_{HHMMSS}.log` (time component is required so
  same-day reruns don't overwrite each other)
- Use `tee` to capture output: `./script.sh 2>&1 | tee logs/script_$(date +%Y%m%d_%H%M%S).log`

## maintenance

- **Audit regularly**: Check for old/orphaned files in `work/`
- **Commit often**: All changes to tracked files should be committed
- **Clean work/**: Remove files that don't follow naming conventions
- **Update readme**: Document new scripts, data sources, or conventions
- **Check .gitignore**: Ensure large/generated files stay untracked

## sealed-data policy

`data/` and any directory marked with `.claude-sealed` are off-limits
to Claude reading content (Claude can still run code that processes
them). For full policy see `~/.claude/CLAUDE.md` and
`~/projects/sealed-data-policy/`.

For metadata sheets containing names or other PII, produce a redacted
version with `redact_metadata.py` and store it in
`work/metadata-redacted/` (gitignored, Claude-readable).
