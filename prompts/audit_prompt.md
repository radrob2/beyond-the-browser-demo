A pipeline in this directory was built by someone else to answer this question: which Hyalella lineages are these samples, and where are their closest known relatives from geographically?

Audit it. You did not build this and should assume nothing in it is correct.

Verify every factual claim in share/report.md against the actual outputs — recompute what can be recomputed and take nothing on trust. Check the pipeline does what docs/conventions.md and the deliverables required. Confirm the QC reference sequence is genuinely absent from the merged dataset, the alignment and the tree — check by sequence content, not only by name or identifier. Check the QC filters do what they claim, including edge cases. Check the nearest-relative table: are the identities and the geography right, and is anything silently dropped or silently missing?

Write TWO things.

1. share/audit.md — the full record. Most serious first, clearly separating what you verified, what you could not verify, and what a domain expert must still check.

2. share/audit.html — a STANDALONE page that will be PROJECTED to a room of 25 people within seconds of you finishing. Do not touch, append to, or modify summary.html. Exactly four screens, each sized to roughly one screenful, in this order:

   Screen 1 — VERDICT. One sentence on whether the report's main biological conclusion survives the audit. Then four counts: claims contradicted, claims misleading or weak, claims verified correct, items needing a domain expert. Then a single pull-quote: your most serious finding, phrased so it can be read aloud verbatim.

   Screen 2 — WHAT'S WRONG. The three most serious findings only. For each: one sentence of claim, one line of evidence with real numbers.

   Screen 3 — WHAT CHECKS OUT. What you recomputed and confirmed. This matters as much as the errors — do not pad it and do not omit it.

   Screen 4 — WHAT NEEDS A HUMAN. What you could not verify, and what only someone who knows these samples can settle.

   Large type, high contrast, self-contained, no external files or network requests, works in light and dark. Assume the back row is five metres away. No paragraph anywhere longer than two sentences.

Open share/audit.html when you are finished.

Do not fix anything. Report only.
