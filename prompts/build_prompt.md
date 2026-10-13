We have cox1 from 30 Hyalella amphipod samples from the Atacama (Sanger PCR), plus cox1 extracted from 10 of our genome assemblies. One Hyalella cox1 reference sequence is for QC only — keep it out of the merged dataset, the alignment and the tree.

The question I want answered: which Hyalella lineages are our samples, and where are their closest known relatives from geographically?

Pull comparable cox1 from public sequence databases — NCBI, BOLD, or both, your call — covering a good range of Hyalella species, plus about 10 other amphipods as outgroups. Keep country and locality metadata wherever it exists; we care about it.

Set up a conda env for this project and install what you need, troubleshooting the env if it fails. Build a working repeatable pipeline following our conventions at docs/conventions.md. Steps: 1. download public seqs. 2. merge with ours. 3. QC: drop seqs under 300 bp; compare to the Hyalella cox1 reference and drop those under 70% identity (reference is for this comparison only); deduplicate identical seqs; ensure unique headers; match our fasta conventions. 4. align with mafft. 5. trim alignment ends, dropping columns where more than 50% of taxa have gaps. 6. build a tree. 7. plot it.

Deliverables, all under share/:

1. The tree as BOTH a PDF and a PNG. It must be readable by someone who has never seen a phylogeny before: colour tip labels by source (public / genome-extracted / pcr), make our own samples visually obvious, and make sure the legend is actually legible — no stray single-letter glyphs as colour keys.

2. A table (TSV) with one row per sample of ours: our sample, its nearest public relative, that relative's species, that relative's country/locality if known, and the % identity.

3. share/report.md — in plain language, for a mixed audience of earth scientists and biologists: what this result means biologically, which of our samples look like known species versus potentially new lineages, where the closest relatives come from, and how confident you are. Full prose belongs here.

4. share/summary.html — a single self-contained page designed to be PROJECTED to a room of 25 people, not read on a laptop. Build it as a sequence of sections, each sized to roughly one screenful, in this order: (a) the question, (b) the tree PNG, as large as fits, (c) a CONDENSED or grouped version of the table that fits one screen at large type — the full table stays in the TSV, (d) five headline findings as short bullets. Large type, high contrast, works in both light and dark. No paragraph longer than two sentences anywhere on this page — link to report.md for the full text. Assume the back row is five metres away.

This is a preliminary v1 and it needs to finish fast. Keep total compute under about five minutes: prefer a fast tree method over exhaustive model selection or heavy bootstrapping, subset where that is defensible, and tell me plainly what you traded away for speed.

Run steps as you complete them so nothing waits until the end. Where work is genuinely independent, dispatch subagents so it runs in parallel. Take advantage of the machine's cores and memory. When everything is finished, open share/summary.html and the tree PDF.
