#!/usr/bin/env python3
"""Step 05: trim alignment ends.

From each end, drop columns in which more than $MAX_GAP_FRACTION of taxa
have a gap, stopping at the first column that passes. Interior columns are
untouched (indels inside the gene are real signal).
"""
import os
import sys

MAX_GAP = float(os.environ["MAX_GAP_FRACTION"])
IN = os.path.join(os.environ["WORK_DIR"], "04_align", "04_align_cox1.fasta")
OUT = os.path.join(os.environ["WORK_DIR"], "05_trim", "05_trim_cox1.fasta")


def main():
    recs, h = [], None
    with open(IN) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                recs.append([line, ""])
            else:
                recs[-1][1] += line.strip()
    ncol = len(recs[0][1])
    n = len(recs)
    gap_frac = [sum(s[i] == "-" for _, s in recs) / n for i in range(ncol)]
    ok = [g <= MAX_GAP for g in gap_frac]
    start = ok.index(True)
    end = ncol - ok[::-1].index(True)
    with open(OUT, "w") as fh:
        for h, s in recs:
            fh.write(f"{h}\n{s[start:end]}\n")
    print(f"[05] {n} sequences; columns {ncol} -> {end - start} "
          f"(trimmed {start} from 5' end, {ncol - end} from 3' end)")


if __name__ == "__main__":
    sys.exit(main())
