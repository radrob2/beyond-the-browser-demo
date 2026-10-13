# Audit of the Atacama *Hyalella* cox1 pipeline and report

*Audited 2026-10-02 against the outputs in `work/`, `share/` and `logs/` from
the run of 2026-10-02 14:13 (with steps 07/09 re-run at 14:23). Nothing in the
pipeline or its outputs was changed. All recomputation was done with separate
scripts writing to a scratch directory, with the same tools (vsearch, IQ-TREE 3,
Biopython) from the `hyalella_dryrun` environment.*

**Bottom line:** the nearest-relative identities in the table are correct, and
so is the broad three-way split (Andean / Titicaca / North American). Several
statements around them are wrong or overstated. The most serious ones: the QC
reference sequence **is** in the dataset, the alignment and the tree. The
support dots on the figure are partly on the wrong branches. The report
over-interprets the "Nevada endemic" and "Titicaca" labels and the
genome-derived results.

---

## A. Wrong: contradicted by the outputs (most serious first)

### A1. The QC reference sequence is in the merged data, the alignment and the tree
The report says (Methods): *"The reference was used only for this filter and is
in no output."* That is false.

- `data/reference/hyalella_cox1_reference.fa` (610 bp) is **byte-identical to
  public record PX696973.1** (*H. azteca*, laboratory culture "Hya-Lyon1",
  Lyon, France, from a toxicity-testing paper). It is also identical to that
  record's BOLD mirror **GBAAX85597-26**.
- Both copies came in through the step-01 download, so the reference sequence
  is present in:

  | file | present as |
  | --- | --- |
  | `work/02_merge/02_merge_cox1.fasta` | PX696973.1 and GBAAX85597-26 |
  | `work/03_qc/03_qc_cox1_all.fasta` | GBAAX85597-26 (PX696973.1 collapsed into it as a duplicate) |
  | `work/03_qc/03_qc_cox1.fasta` (tree set) | GBAAX85597-26 |
  | `work/04_align/04_align_cox1.fasta` | GBAAX85597-26 |
  | `work/05_trim/05_trim_cox1.fasta` | GBAAX85597-26 |
  | `work/06_tree/06_tree_cox1.treefile` | tip GBAAX85597-26 |
  | nearest-relative database (step 03) | both |

- Two further Lyon records (PX696974.1 and PX696975.1, plus their GBAAX
  mirrors) are exact sub-strings of the reference.
- Step 02 only checks that the reference *file* is not concatenated. The header
  `hyalella_cox1_reference` indeed appears nowhere. Nobody checked the sequence
  itself.
- **Consequence:** the nearest public relative of PCR samples **crsp (99.2 %)
  and easy (98.3 %) is the QC reference sequence itself.** If the lab's
  reference or positive control is this Lyon culture, it is a direct candidate
  source of contamination for these two samples. The report does not mention
  this.

### A2. The tree figure puts several support dots on the wrong branches
`07_plot.py` re-roots the IQ-TREE tree with Bio.Phylo `root_with_outgroup`,
twice, and then reads support from each node's label. Re-rooting reverses the
nodes on the path to the new root, but the labels stay on the node objects.
Those labels describe the branch *above* each node, so after re-rooting they
describe the wrong bipartition.

- I compared every node's plotted label with the support of its actual
  bipartition in the IQ-TREE output. **20 of 233 internal nodes carry another
  branch's support.**
- **5 nodes get a black "≥ 95 %" dot although their true UFBoot is below 95.**
  The real values are 45, 77, 83, 59 and 81. The worst is a 199-tip backbone
  node, shown as 99 but truly 45.
- The clade around the genome-derived sequences (61 tips) is plotted as
  supported (90.2/95). Its true support is **5.4/59**.

### A3. "81–95 % bootstrap for the genome-derived sequences and their North American relatives" does not match the tree
- The 10 genome-derived sequences **do not form a group.** Their smallest
  common clade has 61 tips, includes 6 of our PCR samples (crsp, easy, hope,
  pine, urge, xena), and has UFBoot 59.
- For each genome sequence, the smallest clade joining it to a public relative
  has UFBoot 100, 99, 92, 100, 99, 97, 93, 100, **30** and 100 (ata01–ata10).
  I could not find "81–95 %" anywhere. It may have been read off the
  mis-labelled figure (A2).
- The other tree numbers check out: Andean group UFBoot 99 (SH-aLRT 97.4);
  Nevada group 100; Canadian group 99; French group 100.

### A4. "Compared directly with all 1,949 public sequences" uses the wrong number
- 1,949 is the count of distinct sequences **including** our 40 and the 10
  outgroups. Distinct public sequences number **1,899**.
- The nearest-relative search actually used the **3,043** QC-passed public
  sequences *before* deduplication, which is the more complete choice.
- The same wrong number appears again under Limitations: "The
  nearest-relative table always uses all 1,949 sequences".

### A5. "BOLD (2,088 additional records not in GenBank)" and "BOLD records already in GenBank were not duplicated" are false
- Of the 2,088 "BOLD-only" records, **180 are byte-identical to a GenBank
  record that was also downloaded.** 139 have BOLD's `GB…` prefix, meaning BOLD
  mined them from GenBank in the first place.
- The headline "3,109 public records" therefore double-counts some records.
  Exact duplicates are collapsed in QC, so the tree and table are not inflated.
- Deduplication then prefers the **BOLD mirror over the GenBank original**:
  both have a country, so the alphabetical tiebreak wins and `GB…` sorts before
  `PX…`. The mirror has less metadata (no locality, no study title). That is
  why the reference copy in the tree is GBAAX85597-26, not PX696973.1.
- Also note: step 01 reports *"GenBank records given country/locality from
  BOLD: 0"* out of 947 mirrored records. The README's promised metadata
  back-fill never happened. I could not tell whether that is a code problem or
  genuinely no gaps to fill without re-querying BOLD (see C).

### A6. The "Nevada thermal-spring endemic" claim is not supported by the matches
- hope, pine, urge and xena are **exactly as close** to 8 *H. azteca* records
  without any country (DQ464656–59 and their BOLD mirrors GBH058–061, from
  *"DNA barcoding reveals extraordinary cryptic diversity… desert spring
  conservation"*) as they are to the two Hot Creek Spring records (identities
  equal to 0.1 %).
- All of these sit in the same BIN, BOLD:AAA1725 (8 *H. azteca* + 2 "sp. 1").
  The table reports the Nevada record only because of an alignment-length
  tiebreak.
- "*Hyalella* sp. 1 AGC-2024c" is a GenBank placeholder name. That it is the
  "recently described" species, or a "narrow-range spring endemic", is not
  shown by anything in the data.

### A7. Calling *H.* sp. C2 "a Titicaca lineage" is contradicted by the pipeline's own data
- tarn's best hit MN582348.1 (*H.* sp. C2) carries isolate code **3TK10**. That
  is the same population code as the Peruvian "3TK10-1/2/3" records, which sit
  in the *H.* sp. O BIN (ACL3635).
- In the tree, MN582348.1 falls **inside** the 24-tip *H.* sp. O clade.
- "Titicaca" for this record comes only from its publication title, via a
  keyword rule in `08_table.py`.
- isle and silt are also tied (equal identity) between *H.* sp. O and *H.*
  sp. C2.
- The report's treatment of tarn is reasonable. But summary.html and the table
  list "*H.* sp. C2 → Lake Titicaca, Peru/Bolivia" as a separate lineage, which
  is misleading.

### A8. "The three [Titicaca] samples form their own small group in the tree" is not supported
- alba, loch and wand do form a 3-tip clade, but its support is **SH-aLRT 0 /
  UFBoot 55**, so the group itself is unsupported, not only its attachment.
- The "High for the flock" confidence rests entirely on identity (98.1–99.1 %
  to *H. nefrens* / *H.* sp. A1 / Bolivian Titicaca BOLD records), not on the
  tree.

### A9. Smaller factual errors
- "GenBank [Titicaca] records carry no country": OM513737.1 (*H. nefrens*)
  records **Peru**. It is loch's nearest record with a country (98.9 %), so
  "BOLD records from the Bolivian side… 98.1–98.9 %" mixes in a Peruvian
  GenBank record. Bolivian BOLD records are also at 98.9 % for loch, so the
  range survives.
- Genome-derived row, "eastern North America; relatives with locality: Ontario,
  Ohio, **Alberta**": Alberta is western Canada. "Eastern North America" comes
  only from a publication title, and the top hits (JX446358–67) have no
  locality.
- crsp/easy, "also Germany": Swiss records (OP268494.1, 99.1 % / 98.2 %) are
  closer than any German one (98.4 % / 97.4 %). Neither is mentioned in the
  table.
- summary.html: "All 10 genome-derived cox1 are North American ***H.
  azteca***." ata06's nearest is *H.* sp. clade1 JL-2016. The ten span **four
  BINs** (AAA1702 ×7, AAA1706, ACF1948, ACF2658), i.e. several species-level
  lineages of the *H. azteca* complex, not one species.

---

## B. Misleading or weak, though not strictly false

### B1. The 70 % identity QC filter cannot separate *Hyalella* from other amphipods
- **All 10 outgroups pass it**: Caprella 71.1 %, others 73–77 %.
- Our own *Hyalella* samples score 75.0–83.3 % against the reference. Only
  crsp and easy are higher (99.2 % and 98.3 %, see A1).
- So a non-*Hyalella* amphipod mislabelled as *Hyalella* in GenBank or BOLD
  would pass this filter. It only removes non-COI or non-crustacean junk. It
  did that correctly here: all 58 dropped records are only 53–71 % identical
  even to our full-length *Hyalella* cox1, and one is a GenBank record from a
  dung-beetle paper. The report's wording ("not identical enough to our
  *Hyalella* reference") overstates what the filter can do.
- The reference is a single lab-culture *H. azteca* haplotype, which is a
  questionable choice for screening South American taxa.

### B2. QC filter edge cases (verified behaviour)
- **Length (≥ 300 bp):** counts every character, including N and IUPAC codes.
  No passed sequence falls below 300 when only A/C/G/T are counted, so this has
  no effect here. 13 passed sequences have more than 5 % ambiguous bases; the
  maximum is 21.9 %.
- **Identity:** vsearch ignores terminal gaps, and the overlap length is never
  checked. 23 passed sequences overlap the reference by fewer than 300 columns,
  the minimum by 161. Identity is printed to 1 decimal, so 69.95 % would round
  to 70.0 % and pass. The lowest value that passed is 70.1 %, so this had no
  effect here.
- **Reverse-complement:** only sequences that pass the filter are flipped. That
  is correct; 1 was flipped.
- **Header uniqueness:** duplicates get `_2`, `_3` suffixes. A renamed record
  would then fail the accession regex, be classed `unknown`, and lose its
  metadata. This is latent only: 0 renames happened in this run.
- **Dedup:** verified that no identical sequences remain in
  `03_qc_cox1_all.fasta`. None of our 40 sequences is identical to another, or
  to any public record (maximum 99.8 %).
- **End-trimming (step 05):** behaves as documented (1,570 → 638 columns). If no
  column passes, `ok.index(True)` would crash with a bare ValueError. After
  trimming, our samples keep 552–637 bases. Three TWAMI records keep only 294
  of 638 positions (IQ-TREE warns: > 50 % gaps).

### B3. The nearest-relative table silently drops information beyond 25 hits
Step 03 keeps only the top 25 hits per sample. I re-ran the search with up to
5,000 hits:
- **Top hit, identity, aligned length and tie counts all reproduce exactly**
  for all 40 samples.
- **But "closest Chilean record: none among 25 closest" hides real Chilean
  relatives** for alba (95.4 %), loch (96.1 %) and wand (96.0 %). For the North
  American-like samples the closest Chilean record is about 77–80 %; reporting
  that as "none" is fine.
- The "lineage" assigned to each sample is the single top hit, and ties are
  broken by alignment length. Several samples are tied, or within 0.3 %, across
  different names: tarn, isle, silt (sp. O vs sp. C2); alba (13 tied: *H.
  nefrens*, *H.* sp. A1, unnamed Bolivian); hope, pine, urge, xena (sp. 1 vs *H.
  azteca*). The grouped table and the summary page present these as if they
  were clean assignments.
- Identity ignores terminal gaps. kelp's 99.0 % top hit covers only 396 bp,
  and tarn's 98.6 % covers 485 bp. Longer overlaps elsewhere give slightly
  lower identities, which can flip the top hit.
- fern's closest Chilean record is in the Atacama itself (Río Salvador,
  Chacance, **98.5 %**). The report only mentions Chacance for biri (96.5 %).

### B4. Publication-title geography loses its caveat downstream
`08_table.py` correctly labels title-derived locations "(from publication
title)". But `07_plot.py` and `09_summary.py` **strip that qualifier**, so the
figure and the projected page show "eastern North America, lab or wild",
"laboratory culture used in toxicity tests" and "Lake Titicaca, Peru/Bolivia"
as if they were recorded localities. "Laboratory culture" for KU756613.1 is
inferred from the word "toxicity" in a title that does not say "laboratory
culture".

### B5. The contamination interpretation does not fit the genome-derived data well
- The report ranks "contamination from *H. azteca* cultures" as the most likely
  explanation. But the 10 genome-derived sequences are **mutually very
  divergent**: 79.6–99.4 % pairwise. ata08 is about 80 % identical to every
  other one, and ata01 about 92 %. They span four BINs.
- Contamination from a single culture would give one haplotype.
- Instead, each genome sequence matches a *different* record from one GenBank
  study (Major et al., JX4463xx, the lab-vs-wild *H. azteca* paper) at
  99.4–99.8 % over the full 1,560 bp.
- That pattern looks more like the assemblies or extractions being linked to
  that reference panel (sample/file mix-up, or wrong input files) than like
  field or lab contamination. This needs a domain expert and the assembly
  provenance (see D).

### B6. The tree constraint
- **Justification:** the report says that without the constraint "some
  outgroups were pulled inside *Hyalella*". No log in `logs/` documents this.
  I re-ran IQ-TREE without the constraint, and the claim **holds**; see C. The
  constraint costs about 64 log-likelihood units (−32,750.4 vs −32,686.5),
  which the report should state.
- **Circularity:** the constraint puts all 40 of our samples inside *Hyalella*
  by construction. The tree therefore cannot test whether any sample is
  *Hyalella* at all. Identities of 97 % or more make this moot here, but it
  should be stated.
- A deep clade of Canadian *H. azteca* (Manitoba/Ontario BOLD records, plus the
  gappy TWAMI records) sits at the base of *Hyalella* next to the outgroups, so
  *H. azteca* is not monophyletic in this tree.

### B7. Reproducibility and conventions (`docs/conventions.md`)
- `src/`, `config.sh`, `environment.yml`, `readme.md` and all of `share/` are
  **untracked** in git. Conventions: "tracked… commit often".
- `07_plot.py` was edited after the full run (14:23). The figure and
  summary.html come from that later, standalone run, and the
  `run_pipeline` log does not describe them (7 vs 8 bracket groups).
  `report.md` (14:21) predates the final figure.
- Step **08 runs before 07**, which breaks "numbered by pipeline order". It is
  documented, but the numbers no longer mean order.
- `environment.yml` pins no versions, and lists `blast` and `seqkit`, which are
  never used. Results such as the IQ-TREE tie behaviour and the MAFFT output
  may change between versions.
- `share/report.md`, `share/summary.html` and `share/audit.md` do not follow
  the `{step}_{description}` naming pattern (minor).
- The rest checks out: scripts source `config.sh`; bash uses
  `set -euo pipefail`; E-utilities is called over direct HTTP with the API key
  when set; logs are timestamped; `data/` is never modified; directories
  holding our sequences carry `.claude-sealed`.
- The original task brief is not in the repo, so I could only check
  deliverables against `readme.md`. Those are all present: tree PDF/PNG,
  per-sample and grouped table, summary page, report.

---

## C. Verified correct (recomputed)

- Download/QC counts: 1,021 GenBank + 2,088 BOLD-only = 3,109 public
  *Hyalella*; 2,409 with country; 10 outgroups; 30 PCR + 10 genome; 3,159
  merged; 8 dropped < 300 bp (all public); 58 dropped < 70 % (all public); 1
  reverse-complemented; 1,144 duplicates; 1,949 distinct; 236 in the tree.
- All 40 of our sequences pass QC (PCR 552–650 bp, genome 1,560–1,569 bp).
- **Every nearest-relative identity, top hit, aligned length and tie count in
  `08_nearest_relatives.tsv`** reproduces exactly from an independent
  uncapped vsearch search.
- Report numbers: 18 Andean PCR samples (17 sp. O + tarn); sp. O 97.2–99.4 %;
  Menique 97.2–98.7 %; biri→Chacance 96.5 %; tarn 98.6 % (C2) vs 98.3 % (sp. O);
  Titicaca trio 98.1–99.1 %; genome 99.4–99.8 %; crsp/easy 98.3–99.2 %;
  kelp/rove/drft 98.3–99.1 %; hope–xena 97.4–99.0 %; min 97.2 %; 36/40 ≥ 98 %;
  borderline = biri, hope, xena, jolt.
- Calalaste and Quebrada del Diablo are both sp. O localities (Quebrada del
  Diablo is a tied hit, never the reported top hit); BIN ACL3635 contains the
  Menique and Chacance (Chile) records and the Peruvian 3TK10 records.
- Tree: the 18 Andean samples form one clade of 24 tips (with sp. O
  Argentina, Chile, Peru 3TK10-1, and sp. C2), SH-aLRT 97.4 / UFBoot 99. Nevada
  6-tip clade 100/100; Canada 9-tip clade 87.7/99; France 5-tip clade
  (France + Germany) 99.9/100. *Hyalella* + ours vs outgroups 100/100
  (constrained).
- Best model TIM+F+I+R5; 1,570 → 638 columns; log-likelihood −32,750.44.
- BIN AGT3816 (*H. nefrens*) also carries the names *H. cuprea* and *H.
  montforti*.
- **Unconstrained re-run** (same trimmed alignment, TIM+F+I+R5, 1,000 UFBoot,
  seed 12345):
  - *Hyalella* is **not** monophyletic. However the tree is rooted on any one
    outgroup, 7–9 of the other outgroups fall inside *Hyalella*, so the stated
    reason for the constraint is confirmed.
  - All conclusions about our samples are **unchanged** without the
    constraint: Andean 24-tip clade UFBoot 99; Nevada 100; Canada 99; France
    100; Titicaca trio 61 (still unsupported); genome-derived 10 still not a
    group.

---

## D. Not verifiable from the repo: what a domain expert must check

1. **Provenance of the QC reference.** Was `hyalella_cox1_reference.fa` taken
   from PX696973.1, or from a lab positive control or culture? If the lab
   handles the Lyon *H. azteca* culture, crsp and easy (98–99 % to it) are
   prime contamination suspects.
2. **Provenance of the genome assemblies.** Ten divergent sequences, each
   matching a different JX4463xx record over the full gene, need checking
   against the assembly sources and file names before any biology (B5).
3. **Are JX446358/60/61/63/65/67 the lab strain or wild strains** in Major et
   al.? Their GenBank records here carry no strain or locality. The answer
   changes the "lab culture" story for the genome data.
4. **Taxonomy of *H.* sp. C2 / 3TK10 vs *H.* sp. O**, and whether sp. O is now
   described (Puna paper). Also whether *H.* sp. 1 AGC-2024c is the new species
   from the Hot Creek paper or a separate lineage. And where the DQ4646xx
   desert-spring *H. azteca* (same BIN) were collected.
5. **Locality accuracy** of the public records ("Menique lagoon" = Laguna
   Miñiques?; HM138031.1 "Lake Erie, Valens Conservation Area").
6. **BOLD metadata fill = 0.** Re-query BOLD for a few mirrored GenBank
   accessions that lack a country (e.g. MW714609.1) to see whether BOLD has a
   country the code failed to copy.
7. **Outgroup choice.** Check that the 10 outgroup records are correctly
   identified cox1 (several are recent accessions: PZ…, QB…). Also check
   whether Hyalidae/Talitridae/Gammaridae/Caprellidae is a sensible rooting
   for cox1, given saturation.
8. **Biological plausibility.** No two of our 30 PCR samples share a haplotype
   (pairwise 95–98.6 % within the Andean group). Whether that is plausible
   depends on how many sites and individuals per site were sampled, which is
   not recorded in the repo.
9. **Whether 97–98 % cox1 identity means "same species"** in *Hyalella*, where
   cryptic diversity is extreme. The 98 / 96 % bands in `08_table.py` are
   generic, not calibrated for *Hyalella*.
