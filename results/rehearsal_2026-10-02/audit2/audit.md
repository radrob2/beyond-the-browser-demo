# Audit: Atacama *Hyalella* cox1 pipeline and `share/report.md`

*Audited 2026-10-02 against the outputs in `work/` and `share/` from the 14:13–14:23 run (logs `logs/*_20261002_14*.log`).
Nothing was fixed or re-run in the project. All recomputation used separate scripts and scratch copies outside the project.*

## Verdict

**The report's main biological conclusion holds up.** 18 PCR samples belong to the unnamed high-Andean
*H.* sp. O lineage (BOLD:ACL3635: Argentina, Chile, Peru). 3 PCR samples sit with the Titicaca flock. 19
sequences (all 10 genome-derived and 9 PCR) match North American *H. azteca*-complex lineages, and no sample is
below 97.2 % identity. I reproduced every nearest-relative identity exactly.

**Several of its supporting statements are wrong, and one of them bears directly on the contamination question.**
The QC reference sequence, which the report says "is in no output", is in the alignment, the tree and the figure.
It is also the nearest public relative of PCR samples **crsp** and **easy**.

| | count |
| --- | --- |
| Report claims contradicted | **6** |
| Report claims misleading or weak | **6** |
| Report claims verified correct | **30** |
| Items needing a domain expert | **6** |

The per-claim ledger is in the "Claim-by-claim ledger" section below.

---

## A. Findings, most serious first

### 1. The QC reference is in the merged data, alignment, tree and figure, and it is the nearest relative of two samples  [CONTRADICTED]

*Report, Methods: "The reference was used only for this filter and is in no output."*

I compared sequence content, not names: exact match, substring in both orientations, and 40-mer sharing.

- `data/reference/hyalella_cox1_reference.fa` holds one 610 bp sequence. It is **identical over all 610 bp** to
  GenBank **PX696973.1** (*H. azteca*, isolate "Hya-Lyon1", a Lyon laboratory culture from a toxicity-test study).
  It is also identical to that record's BOLD copy, **GBAAX85597-26**, whose sampleid field reads "PX696973".
- Where the reference sequence appears (exact 610/610 match):

  | file | present as |
  | --- | --- |
  | `work/01_download/*_ncbi.fasta` | PX696973.1 |
  | `work/01_download/*_bold.fasta` | GBAAX85597-26 |
  | `work/02_merge/02_merge_cox1.fasta` | both |
  | `work/03_qc/03_qc_cox1_all.fasta` | GBAAX85597-26 (PX696973.1 collapsed into it as a duplicate) |
  | `work/03_qc/03_qc_cox1.fasta` (tree set) | GBAAX85597-26 |
  | `work/04_align`, `work/05_trim` | GBAAX85597-26 |
  | `work/06_tree/06_tree_cox1.treefile` | tip GBAAX85597-26 |
  | `share/figure/07_tree_cox1.pdf` | labelled tip "Hyalella azteca \| France \| GBAAX85597-26" |
  | `share/08_nearest_relatives.tsv` | PX696973.1 = nearest relative of **crsp** (99.2 %, 610 bp) and **easy** (98.3 %, 574 bp) |

- Two sister lab-culture haplotypes (PX696974.1 and PX696975.1, plus their BOLD copies GBAAX87964-26 and
  GBAAX87589-26) are exact substrings of the reference. GBAAX87589-26 is also in the tree.
- **Why it happened:** `02_merge.sh` leaves out the reference *file* by path. Nothing checks whether the same
  sequence arrives through the public download.
- In the tree, crsp sits in a 5-tip clade with GBAAX85597-26, GBAAX87589-26 and a German record. That clade has
  SH-aLRT 99.9 and UFBoot 100.
- **Why this matters:** the report's top hypothesis for the North American matches is contamination. The lab's own
  QC reference is a lab-culture *H. azteca* sequence, and two PCR samples match it best. Whether this lab handles
  that culture or its DNA is a question for the people who ran the samples (see section C).

### 2. The tree figure puts support values on the wrong nodes, and the report quotes the wrong values  [CONTRADICTED]

*Report: "81–95 % for the genome-derived sequences and their North American relatives". Report and figure: black dots = UFBoot ≥ 95.*

- IQ-TREE wrote the tree rooted at tip ata01. `07_plot.py` re-roots it with Biopython's `root_with_outgroup`, which
  leaves each support label on its original node. On every node along the old-root-to-new-root path, the label then
  describes the wrong bipartition.
- I recomputed support from bipartitions of the unrooted treefile. **20 of 233 internal nodes** carry the wrong
  label after re-rooting. On **10** of them the "strong support" dot is wrong (shown where support is < 95, or
  missing where it is ≥ 95).
- The report's "81–95 %" are exactly the shifted labels:

  | clade | plotted (SH-aLRT/UFBoot) | true |
  | --- | --- | --- |
  | 8 genome samples + relatives (23 tips) | 99.8/81 | **69.3/74** |
  | all 10 genome samples (61 tips; also contains crsp, easy, hope, pine, urge, xena) | 90.2/95 | **5.4/59** |
  | 9 genome samples + relatives (26 tips) | 69.3/74 | **99.8/98** |
  | ata01 + JX446365.1 (2 tips) | none | **88.6/100** |

- True UFBoot for the smallest clade joining each genome sample to a public sequence ranges from **30 (ata09) to
  100**. It is not 81–95.
- Not affected: the sp. O clade (97.4/99), Nevada (100/100), Alberta (87.7/99), France (99.9/100) and Titicaca
  (0/55). Those nodes are not on the re-rooting path, and the report's values for them are correct.

### 3. "BOLD records already in GenBank were not duplicated" is false  [CONTRADICTED]

*Report: "BOLD, … (2,088 additional records not in GenBank)". Readme and Methods say the same.*

- **139** of the 2,088 "BOLD-only" records have BOLD's `GB…` process IDs, meaning BOLD mined them from GenBank.
  **132** of these are sequence-identical to a GenBank record in the same download. Overall, 180 BOLD-only records
  are identical to some downloaded GenBank sequence.
- Cause: `01_download.py` detects mirrors only through BOLD's `insdc_acs` field, which is empty for these
  records. GBAAX85597-26, for example, carries the accession only in `sampleid`.
- Consequences:
  - "3,109 public records" counts some specimens twice.
  - `n_equally_close` is inflated. crsp's "4 equally close" are 2 sequences × 2 copies, and the Nevada samples'
    "4" are 2 + 2 (DQ464656.1 = GBH058-06).
  - In the tree, a geography-less BOLD copy is shown instead of the GenBank record that has a locality.
    Deduplication breaks ties alphabetically, so GBAAX85597-26 beats PX696973.1 "Lyon".
- Within-pipeline deduplication by identical sequence removes most of the damage from the tree. Counts and
  tie counts remain affected.

### 4. The Nevada "thermal-spring species" relative is a tie-break, not a unique nearest match  [MISLEADING]

*Report: hope, pine, urge, xena → "Hyalella sp. 1 AGC-2024c, Hot Creek Spring, Nye County, Nevada (a recently described thermal-spring species)"; "a narrow-range spring endemic".*

- For all four samples, **DQ464656.1 (*H. azteca*, no country or locality in the record)** is exactly as close as
  the Hot Creek records (97.4 / 98.6 / 99.0 / 97.7 %, all over 621 bp). DQ464656.1 comes from Witt et al.'s
  desert-spring barcoding study.
- The table shows the Nevada record only because sorting among exact ties falls that way.
- BIN BOLD:AAA1725 holds 10 records. 8 of them are labelled *H. azteca* and have no country.
- The GenBank organism name is "*Hyalella* sp. 1 AGC-2024c", an informal label. A described species or a "narrow
  range" cannot be read from these data. The identity values themselves are correct.

### 5. The genome-derived samples' "eastern North America" comes from a paper title, and their hits are not in the lab-strain BIN  [MISLEADING]

*Report table: all 10 genome-derived → "H. azteca complex | eastern North America; relatives with locality: Ontario, Ohio, Alberta | 99.4–99.8".*

- 9 of the 10 top hits are JX4463xx records, and **none of them records a country**. "Eastern North America" is
  inferred from the paper title ("The common ecotoxicology laboratory strain of *H. azteca* is genetically distinct
  from most wild strains sampled in eastern North America"). The table itself says "lab or wild". The report drops
  that ambiguity.
- The genome hits fall in BINs AAA1702, AAA1706, ACF1948 and ACF2658. The same paper's JX446307–318 records sit in
  **BOLD:AAX2514**, the BIN that also holds this run's "(lab culture)" localities and the Lyon culture. On BIN
  evidence, the genome sequences look like *wild* North American haplotypes, not the standard lab strain. Only the
  PCR samples crsp and easy match the lab-culture BIN.
- Alberta is not eastern North America (ata04, ata05 and ata10 have their nearest located relative in Alberta).
- ata06's top hit is "*Hyalella* sp. clade1 JL-2016" (from a Cu/Ni toxicity study), not *H. azteca* by name.
  The 99.4 % range still holds.

### 6. tarn's "Titicaca lineage" hit is probably the same population as sp. O  [MISLEADING]

*Report: tarn "technically has its single best hit to Hyalella sp. C2 (a Titicaca lineage, 98.6 %)".*

- MN582348.1 (sp. C2) has isolate code **3TK10**. That is the same code as the three Peruvian "*Hyalella* sp.
  3TK10-1/-2/-3" records that sit **inside** BIN ACL3635 (sp. O).
- The ranking depends on the method. With BLAST, sp. C2 is also the top hit for **isle** (98.2 vs 98.0) and
  **jolt** (98.0 vs 97.7). With vsearch, tarn's C2 hit covers only 485 bp, against 589 bp for sp. O.
- The grouped table and `summary.html` still list "Hyalella sp. C2 – Lake Titicaca, Peru/Bolivia" as a separate
  lineage row for tarn. The report's conclusion (tarn belongs with sp. O) is right. The table row is misleading.

### 7. "Compared with all 1,949 public sequences" uses the wrong number  [CONTRADICTED]

- The 1,949 unique QC-passed sequences include our 40 and 10 outgroups, so 1,899 are public.
- The nearest-relative search actually ran against **3,043** QC-passed public records *before* deduplication
  (`03_qc.py` step 5). That is the better design, but the report states it wrongly twice ("all 1,949 public
  sequences", "the nearest-relative table always uses all 1,949 sequences").

### 8. Río Salvador is the closest Chilean record for two samples, not one  [CONTRADICTED]

*Report: "One sample (biri) is closest among Chilean records to an animal from the Río Salvador at Chacance (96.5 %)".*

- **fern** is also closest among Chilean records to TTKK287-11, Río Salvador (98.5 %).
- The 97.2–98.7 % range for the other 16 sp. O samples against Menique-lagoon TTKK299-11 is correct.

### 9. The borderline samples contradict the report's own definition  [MISLEADING]

- The report says within-species cox1 differs by "0–2 %". It then calls biri (97.2 %, i.e. 2.8 % divergent),
  hope (2.6 %), xena (2.3 %) and jolt (2.1 %) "still within the range usually seen inside one species".
- All four lie in the report's own 2–3 % grey zone. The "no new species" conclusion is weaker than stated for
  biri and hope.

### 10. The QC identity filter is weaker than described  [WEAK]

*Report: 58 records "< 70 % identical to our Hyalella reference (most likely mis-labelled or a different gene region)".*

- **All 10 non-*Hyalella* outgroups pass the 70 % filter (71.1–76.6 %).** So the filter cannot catch mislabelled
  amphipods. It only checks "is this cox1".
- In vsearch, `N` counts as a **match**. A synthetic test confirmed that 20 % N scores 100 % identity. 217 merged
  records contain N, and 13 QC-passed public records are > 5 % N (up to 124 N in 566 bp, NJCGS1170-11). None of
  them is a reported nearest relative (all 17 top-hit targets have 0 N), so the results are unaffected.
- Identity ignores terminal gaps, so partial overlaps pass. 15 EU6217xx records align only 161–264 bp to the
  reference and still pass. In the trimmed alignment, every tree sequence has ≥ 294 of 638 sites. IQ-TREE flags 3
  sequences > 50 % gaps (TWAMI095/098/100-09).
- Of the 58 low-identity drops, 39 did not align at all at ≥ 50 % identity. 57 are BOLD records, and **56 are
  Canadian "*H. azteca*"** with no BIN. One GenBank record comes from a dung-beetle paper. "Mis-labelled or a
  different gene region" is plausible but unverified (see section C).
- Length filter: correct. Exactly 8 records < 300 bp were removed, all public BOLD (207–276 bp). No QC-passed
  record has < 300 informative (non-N) bases.

### 11. The nearest-relative table hides some relevant context

- alba, loch and wand: "closest Chilean record: none among 25 closest". The full search finds Chilean records at
  **95.4–96.1 %** (TTKK289-11, Loa River at Chiu Chiu, BIN AHQ7179), just outside the 25-hit window. In an
  Atacama study, that is worth stating.
- crsp and easy: `relative_bin` is blank because the GenBank record has no BIN. Its BOLD copy is in
  **BOLD:AAX2514**, a lab-culture-rich BIN with records from Germany, USA, Mexico, Switzerland and France.
- 14 of 40 top hits show `relative_country = not recorded`. That is handled transparently through the
  "nearest with country" columns.

### 12. Code edge cases that this run did not trigger (latent bugs)

- `08_table.py` lists *every* pcr/genome sample, including ones that failed QC. Such a sample would show
  "no public hit ≥ 70 %" instead of "failed QC". All 40 passed this time.
- `03_qc.py` renames duplicate headers to `NAME_2`. That breaks the accession regex, so the source becomes
  "unknown" and the country lookup fails. No duplicates occurred (0 renamed).
- `03_qc.py` filters on raw length, including N. No record was affected.
- The top-hit choice among exact ties depends on sort order (identity, then aligned length). This is the root of
  finding 4.

### 13. Conventions (`docs/conventions.md`) and deliverables

Checks that passed:

- Sealed markers are present on `work/02_merge` to `work/05_trim`. Nothing writes to `data/`.
- `work/06_tree` (unsealed) contains no sequence text. Its only ≥ 60-character nucleotide-alphabet runs are dash
  separators in the log.
- Bash scripts use `set -euo pipefail` and source `config.sh`.
- Logs follow `{script}_{YYYYMMDD}_{HHMMSS}.log`.
- FASTA outputs are single-line and uppercase (checked steps 02–05).
- Headers are trimmed only to the leading token.
- Downloads use direct E-utilities HTTP with history server and API key.
- Files are named in lowercase.
- The readme is updated.

Deviations:

- The readme's description of step 01 ("BOLD records already in GenBank are not duplicated") is false (finding 3).
- `src/`, `config.sh`, `readme.md` and all of `share/` are **untracked in git**. The conventions say these are
  tracked and should be committed often. No committed version of the code that produced these outputs exists.
- `07_plot.py`, `08_table.py` and `09_summary.py` lack the execute bit. They are called via `python`, so this is
  harmless.

Deliverables:

- The question (lineage and geography of nearest relatives) is answered by `08_nearest_relatives*.tsv`, the tree
  figure and the report.
- The figure carries finding 2. The tables carry findings 3, 4 and 11.
- `summary.html` repeats the "eastern North America" wording of finding 5 and the separate sp. C2 row of finding 6.

---

## B. What I verified (recomputed, not trusted)

| item | method | result |
| --- | --- | --- |
| Download counts: 1,021 GenBank, 2,088 BOLD-only, 3,109 total, 2,409 with country (22.5 % without), 10 outgroups (one per genus) | recount from `01_download_metadata.tsv` | match (but see finding 3 on duplicates) |
| QC: 8 short, 58 low-identity (all public), 40/40 ours pass, 1 reverse-complemented (FJ375302.1), 1,949 unique, 1,144 duplicates | recount from `03_qc_report.tsv` | match |
| Tree set 236 = 186 public + 30 PCR + 10 genome + 10 outgroups; 141 clusters at 95 % | recount | match; tree tips = tree set exactly |
| Alignment 1,570 → 638 columns; ends-only trimming at > 50 % gaps | recount + code read | match |
| IQ-TREE: TIM+F+I+R5 (BIC), 1000 UFBoot + 1000 SH-aLRT, constraint applied, *Hyalella* + ours monophyletic, rooted on outgroups | `.iqtree`, log, treefile | match |
| Nearest relative, identity and tie count for **all 40 samples** | independent untruncated vsearch search of all 3,043 QC-passed public records | **40/40 identical** to the table; the 25-hit cap changed nothing |
| Same, by a different method | BLASTn (HSP ≥ 300 bp) | same organism or BIN for 37/40; same lineage for 38/40 (isle and jolt flip to sp. C2, see finding 6); identities within 0.7 % |
| No QC-dropped record would have been a closer relative | ours vs all 66 dropped records | 0 of 120 hits beat the reported nearest |
| Identity ranges: sp. O 97.2–99.4; Titicaca 98.1–99.1; Alberta 98.3–99.1; Nevada 97.4–99.0; France 98.3–99.2; genome 99.4–99.8; all ≥ 97.2; 36/40 ≥ 98; borderline biri 97.2, hope 97.4, xena 97.7, jolt 97.9 | table + recomputation | match |
| Calalaste and Quebrada del Diablo as closest sp. O localities | ties in full search | Quebrada del Diablo (OR726347.1) tied-closest for 8 samples; Calalaste for 10; 7 samples' only top hit (OR726343.1) has country but no locality |
| BIN ACL3635 spans Argentina (4), Chile (6, Antofagasta incl. "Menique lagoon"), Peru (3, "3TK10") | metadata | match |
| Chilean identities 97.2–98.7 % (Menique), biri 96.5 % (Río Salvador) | full search | match (but fern also Río Salvador, finding 8) |
| Titicaca: Bolivian BOLD 98.1–98.9 %; *H. nefrens* BIN AGT3816 also holds *H. cuprea* and *H. montforti* | full search + metadata | match |
| Tree: 18 sp. O samples form one clade with Argentine, Chilean and Peruvian public records, SH-aLRT 97.4 / UFBoot 99; tarn inside | bipartition support | match |
| Tree: Nevada 100, Alberta 99, France 100 UFBoot; Titicaca trio separate, attachment UFBoot 55 | bipartition support | match |
| Genome cox1 lengths ~1,560 (1,560–1,569); PCR 552–650 bp; ≤ 4 ambiguous bases per sample | length/composition counts | match |
| No sequence text in unsealed `work/06_tree/` | pattern scan | clean |

## C. What a domain expert must still check

1. **Does this lab use the Lyon *H. azteca* culture sequence (Hya-Lyon1, PX696973.1) as a reference or a positive
   control, or handle that culture?** crsp (99.2 %) and easy (98.3 %) match it best. This is the most direct
   contamination lead in the data.
2. **Genome assemblies:** were the 10 genome animals collected at the same Atacama sites as the PCR animals?
   And are the JX4463xx hits (Major et al.) lab or wild animals? On BIN evidence they look wild (finding 5).
3. **Nevada matches:** do the tied *H. azteca* records from the desert-spring study (DQ464656–59) come from the
   same Great Basin springs as Hot Creek? Is "*Hyalella* sp. 1 AGC-2024c" now a described, narrow-range species?
4. ***H.* sp. C2 / "3TK10":** is this a Titicaca-flock lineage, or the same population as sp. O? This decides
   whether tarn (and arguably isle and jolt) need any Titicaca caveat.
5. **The 56 Canadian BOLD "*H. azteca*" records** that fail the 70 % filter: NUMTs, a non-barcode region,
   contamination or mislabelling? They do not affect our nearest relatives, but the report's explanation is a guess.
6. **Taxonomy:** whether sp. O, at 2–3 % from biri and jolt, warrants description. This needs morphology and
   additional loci.

## D. What I could not verify

- **"Unconstrained, some outgroups were pulled inside *Hyalella*."** No unconstrained tree was kept. The only
  earlier run whose constraint failed ("vs 0 outgroups", log `06_tree_20261002_140505.log`) was overwritten.
  Checking this needs a re-run, which I did not do.
- **Completeness of the downloads.** I did not re-query NCBI or BOLD, and the raw BOLD TSV was not retained. So I
  cannot tell why 0 GenBank records received country data from BOLD, even though 254 BOLD-mirrored GenBank records
  lack a country.
- **Publication-title geography** ("Lake Titicaca, Peru/Bolivia", "eastern North America"). These are inferences
  from paper titles, not record metadata, and I did not read the papers.
- **MAFFT and IQ-TREE reproducibility.** I did not re-run them. I checked their outputs for internal consistency
  only.
- **The PNG figure and `summary.html` visuals.** I checked the PDF's text layer and the text of `summary.html`, not
  the rendered images.

---

## Claim-by-claim ledger (`share/report.md`)

| # | claim | verdict |
| --- | --- | --- |
| 1 | 40 sequences fall into three groups; 18/30 PCR high-Andean | verified |
| 2 | That lineage known from Argentine Puna, Antofagasta (Chile), Peru | verified |
| 3 | 3 PCR samples sit with the Titicaca flock | verified |
| 4 | All 10 genome + 9 PCR almost identical to North American lineages | verified |
| 5 | …including lineages kept as laboratory cultures | verified (crsp/easy) |
| 6 | Every sample ≥ 97 % to a public sequence | verified |
| 7 | GenBank 1,021 records | verified |
| 8 | BOLD 2,088 records "not in GenBank"; BOLD records in GenBank not duplicated | **contradicted** (F3) |
| 9 | 3,109 public, 2,409 with country | verified (counts) |
| 10 | 10 outgroup genera | verified |
| 11 | QC removed 8 < 300 bp and 58 < 70 % | verified |
| 12 | The 58 are "most likely mis-labelled or a different gene region" | **weak** (F10) |
| 13 | All 40 of ours passed; 1,949 distinct sequences | verified |
| 14 | Compared with "all 1,949 public sequences" | **contradicted** (F7) |
| 15 | Tree subset design, 236 sequences | verified |
| 16 | 18 samples one tight group, bootstrap 99 %, with AR/CL/PE records | verified |
| 17 | Closest relatives sp. O (ACL3635), Calalaste / Quebrada del Diablo, 97.2–99.4 % | verified |
| 18 | Cluster includes Antofagasta ("Menique lagoon"); 97.2–98.7 % to Chilean animals | verified |
| 19 | Only biri is closest to Río Salvador | **contradicted** (F8) |
| 20 | Peruvian "3TK10" specimens in the cluster | verified |
| 21 | tarn: best hit sp. C2, "a Titicaca lineage" | **misleading** (F6) |
| 22 | alba, loch, wand 98.1–99.1 % to *H. nefrens* / sp. A1 | verified |
| 23 | Bolivian BOLD records 98.1–98.9 % | verified |
| 24 | Titicaca trio own group; attachment poorly supported | verified |
| 25 | *H. nefrens* BIN also holds *H. cuprea*, *H. montforti* | verified |
| 26 | Genome-derived: "eastern North America" | **misleading** (F5) |
| 27 | crsp, easy: Lyon lab culture, also Germany, 98.3–99.2 % (omits that this *is* the QC reference) | **misleading** (F1) |
| 28 | kelp, rove, drft: Alberta, 98.3–99.1 % | verified |
| 29 | hope/pine/urge/xena: Nevada thermal-spring species, narrow endemic | **misleading** (F4) |
| 30 | Nevada, Canadian, French groups 99–100 % bootstrap | verified |
| 31 | Genome-derived groups 81–95 % bootstrap | **contradicted** (F2) |
| 32 | 36 of 40 ≥ 98 %; 4 borderline with stated values | verified |
| 33 | Borderline samples "within the range usually seen inside one species" | **misleading** (F9) |
| 34 | ~640 bases / 638 aligned positions; genome ~1,560 bases | verified |
| 35 | About a quarter of public records lack a country | verified (22.5 %) |
| 36 | Unconstrained tree pulled outgroups into *Hyalella* | could not verify |
| 37 | One public sequence reverse-complemented | verified |
| 38 | Reference used only for the filter and "in no output" | **contradicted** (F1) |
| 39 | Exact duplicates removed, ours kept first | verified |
| 40 | Trim > 50 % gaps 1,570 → 638; IQ-TREE TIM+F+I+R5, UFBoot + SH-aLRT, constraint, outgroup rooting | verified |
| 41 | vsearch exhaustive global search; identity = matches / aligned columns, terminal gaps ignored | verified (N counts as match, F10) |
| 42 | Figure: black dots mark UFBoot ≥ 95 | **contradicted** (F2: 10 dots wrong) |
| 43 | Methods: data from E-utilities (direct HTTP) and BOLD v5 API, downloaded 2026-10-02 | verified (code + logs) |

Totals: **30 verified**, **6 contradicted** (8, 14, 19, 31, 38, 42), **6 misleading or weak** (12, 21, 26, 27, 29, 33),
1 could not verify (36). **6** domain-expert items (section C).
