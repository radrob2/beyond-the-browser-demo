#!/usr/bin/env python3
"""Step 03: QC the merged cox1 set and pick the public subset for the tree.

QC, in order:
  1. drop sequences shorter than $MIN_LEN bp
  2. compare each sequence to the Hyalella cox1 reference (vsearch global
     alignment, both strands; identity = matches / aligned columns, terminal
     gaps ignored); reverse-complement any sequence that hits on the minus
     strand; drop those under $MIN_IDENTITY %
  3. trim headers to the leading identifier token and make them unique
  4. deduplicate identical sequences (keep: genome > pcr > outgroup >
     public with country > public), recording what each kept tip stands for

Then:
  5. nearest public relatives: each of our QC-passed sequences (before
     dedup) vs every QC-passed public Hyalella sequence, exhaustive search
  6. tree set: cluster the deduplicated public Hyalella at $CLUSTER_ID identity
     (centroids preferably records with country and species name); keep
     each centroid plus one extra member per additional South American
     country; always keep the top public hits of each of our samples, all
     outgroups and all of our samples

The reference itself never enters any output FASTA.
"""
import collections
import os
import re
import subprocess
import sys

import pandas as pd
from Bio.Seq import Seq

WORK = os.environ["WORK_DIR"]
OUTDIR = os.path.join(WORK, "03_qc")
MERGED = os.path.join(WORK, "02_merge", "02_merge_cox1.fasta")
META = os.path.join(WORK, "01_download", "01_download_metadata.tsv")
REF = os.environ["REFERENCE_FA"]
MIN_LEN = int(os.environ["MIN_LEN"])
MIN_ID = float(os.environ["MIN_IDENTITY"])
THREADS = os.environ["THREADS"]
TOP_HITS = int(os.environ["TOP_HITS_PER_SAMPLE"])
CLUSTER_ID = float(os.environ["CLUSTER_ID"])

ACCESSION = re.compile(r"^[A-Z]{1,6}_?\d+\.\d+$")   # GenBank, e.g. KU870690.1
PROCESSID = re.compile(r"^[A-Z0-9]+-\d{2}$")         # BOLD, e.g. HYAL123-19
SOUTH_AMERICA = {"Argentina", "Bolivia", "Brazil", "Chile", "Colombia", "Ecuador",
                 "Paraguay", "Peru", "Uruguay", "Venezuela", "Guyana", "Suriname",
                 "French Guiana", "Falkland Islands (Islas Malvinas)"}


def read_fasta(path):
    recs, h = [], None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                h = line[1:]
                recs.append([h, ""])
            elif h is not None:
                recs[-1][1] += line.strip()
    return recs


def write_fasta(path, items):
    with open(path, "w") as fh:
        for h, s in items:
            fh.write(f">{h}\n{s}\n")


def source_of(name, outgroup_ids):
    if name.startswith("genome_"):
        return "genome"
    if "_pcr_" in name:
        return "pcr"
    if name in outgroup_ids:
        return "outgroup"
    if ACCESSION.match(name) or PROCESSID.match(name):
        return "public"
    return "unknown"


def vsearch(*args):
    subprocess.run(["vsearch", "--threads", THREADS, "--quiet", *args], check=True)


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    open(os.path.join(OUTDIR, ".claude-sealed"), "a").close()
    meta = pd.read_csv(META, sep="\t", keep_default_na=False)
    outgroup_ids = set(meta.loc[meta.group == "outgroup", "accession"])
    country = dict(zip(meta.accession, meta.country))

    recs = read_fasta(MERGED)
    print(f"[03] input records: {len(recs)}")

    # 3. header trimming + uniqueness (done first so every table keys on it)
    seen = collections.Counter()
    rows, seqs = [], {}
    for h, s in recs:
        name = h.split()[0]
        seen[name] += 1
        if seen[name] > 1:
            name = f"{name}_{seen[name]}"
        seqs[name] = s
        rows.append({"name": name, "source": source_of(name, outgroup_ids),
                     "length": len(s)})
    qc = pd.DataFrame(rows).set_index("name")
    n_renamed = sum(c - 1 for c in seen.values())
    print(f"[03] headers trimmed to identifier; duplicates renamed: {n_renamed}")
    if (qc.source == "unknown").any():
        print(f"[03] WARNING unknown source: {list(qc.index[qc.source == 'unknown'])}")

    # 1. length
    qc["pass_length"] = qc.length >= MIN_LEN
    print(f"[03] dropped < {MIN_LEN} bp: {(~qc.pass_length).sum()}")

    # 2. identity to reference
    tmp_in = os.path.join(OUTDIR, "03_qc_tmp_lenpass.fasta")
    write_fasta(tmp_in, [(n, seqs[n]) for n in qc.index[qc.pass_length]])
    ref_hits = os.path.join(OUTDIR, "03_qc_reference_identity.tsv")
    vsearch("--usearch_global", tmp_in, "--db", REF, "--id", "0.5",
            "--strand", "both", "--maxaccepts", "1", "--maxrejects", "0",
            "--output_no_hits", "--userout", ref_hits,
            "--userfields", "query+id+qstrand+alnlen")
    hits = pd.read_csv(ref_hits, sep="\t", header=None,
                       names=["name", "ref_identity", "strand", "ref_alnlen"])
    hits = hits.drop_duplicates("name").set_index("name")
    qc = qc.join(hits)
    qc["ref_identity"] = qc.ref_identity.fillna(0.0)
    qc["pass_identity"] = qc.pass_length & (qc.ref_identity >= MIN_ID)
    flipped = qc.index[qc.pass_identity & (qc.strand == "-")]
    for n in flipped:
        seqs[n] = str(Seq(seqs[n]).reverse_complement())
    print(f"[03] reverse-complemented to match reference: {len(flipped)}")
    print(f"[03] dropped < {MIN_ID:g}% identity to reference: "
          f"{(qc.pass_length & ~qc.pass_identity).sum()}")
    for src, g in qc[qc.pass_length & ~qc.pass_identity].groupby("source"):
        print(f"       {src}: {len(g)}")
    os.remove(tmp_in)

    # 4. deduplicate identical sequences
    rank = {"genome": 0, "pcr": 1, "outgroup": 2, "public": 3, "unknown": 4}
    passed = qc[qc.pass_identity].copy()
    passed["priority"] = [(rank[r.source], 0 if country.get(n) else 1, n)
                          for n, r in passed.iterrows()]
    groups = collections.defaultdict(list)
    for n in passed.sort_values("priority").index:
        groups[seqs[n]].append(n)
    qc["kept_as"] = ""
    for members in groups.values():
        for m in members:
            qc.loc[m, "kept_as"] = members[0]
    keep = [m[0] for m in groups.values()]
    print(f"[03] identical-sequence duplicates removed: {len(passed) - len(keep)}")
    qc["status"] = "kept"
    qc.loc[~qc.pass_length, "status"] = "short"
    qc.loc[qc.pass_length & ~qc.pass_identity, "status"] = "low_identity"
    qc.loc[qc.pass_identity & (qc.kept_as != qc.index), "status"] = "duplicate"
    qc.drop(columns=["pass_length", "pass_identity"]).to_csv(
        os.path.join(OUTDIR, "03_qc_report.tsv"), sep="\t")
    write_fasta(os.path.join(OUTDIR, "03_qc_cox1_all.fasta"),
                [(n, seqs[n]) for n in keep])
    print(f"[03] QC-passed unique sequences: {len(keep)} "
          + str(qc.loc[keep].source.value_counts().to_dict()))

    # 5. nearest public relatives of each of our sequences (pre-dedup)
    ours = list(passed.index[passed.source.isin(["pcr", "genome"])])
    public = list(passed.index[passed.source == "public"])
    q_fa = os.path.join(OUTDIR, "03_qc_tmp_ours.fasta")
    db_fa = os.path.join(OUTDIR, "03_qc_tmp_public.fasta")
    write_fasta(q_fa, [(n, seqs[n]) for n in ours])
    write_fasta(db_fa, [(n, seqs[n]) for n in public])
    nn = os.path.join(OUTDIR, "03_qc_nearest_public.tsv")
    vsearch("--usearch_global", q_fa, "--db", db_fa, "--id", "0.7",
            "--strand", "plus", "--maxaccepts", "0", "--maxrejects", "0",
            "--maxhits", "25", "--output_no_hits", "--userout", nn,
            "--userfields", "query+target+id+alnlen+mism+opens")
    os.remove(q_fa)
    os.remove(db_fa)
    nnd = pd.read_csv(nn, sep="\t", header=None,
                      names=["query", "target", "identity", "alnlen", "mism", "gap_opens"])
    nnd = nnd.sort_values(["query", "identity", "alnlen"], ascending=[True, False, False])
    nnd.to_csv(nn, sep="\t", index=False)
    print(f"[03] nearest-relative search: {nnd['query'].nunique()} of ours, "
          f"{len(nnd)} hits")

    # 6. tree set
    # centroid preference: has country > named species > longer
    organism = dict(zip(meta.accession, meta.organism))
    named = lambda n: " sp." not in organism.get(n, "") and organism.get(n, "") != "Hyalella"
    pub_keep = sorted((n for n in keep if qc.loc[n, "source"] == "public"),
                      key=lambda n: (not country.get(n), not named(n), -len(seqs[n]), n))
    cl_in = os.path.join(OUTDIR, "03_qc_tmp_cluster.fasta")
    uc = os.path.join(OUTDIR, "03_qc_public_clusters.uc")
    write_fasta(cl_in, [(n, seqs[n]) for n in pub_keep])
    vsearch("--cluster_smallmem", cl_in, "--usersort", "--id", str(CLUSTER_ID),
            "--strand", "plus", "--uc", uc)
    os.remove(cl_in)
    clusters = collections.defaultdict(list)
    with open(uc) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == "S":
                clusters[f[8]].insert(0, f[8])
            elif f[0] == "H":
                clusters[f[9]].append(f[8])
    selected = set()
    for centroid, members in clusters.items():
        selected.add(centroid)
        countries = {country.get(centroid, "")}
        for m in members[1:]:
            c = country.get(m, "")
            if c in SOUTH_AMERICA and c not in countries:
                selected.add(m)
                countries.add(c)
    print(f"[03] public {CLUSTER_ID:.0%} clusters: {len(clusters)}; representatives (centroid "
          f"+ 1 per extra South American country): {len(selected)}")
    n_before = len(selected)
    for q, g in nnd[nnd.target != "*"].groupby("query"):
        reps = list(dict.fromkeys(qc.loc[g.target, "kept_as"]))  # in identity order
        selected.update(reps[:TOP_HITS])
    print(f"[03] added as nearest relatives of our samples: {len(selected) - n_before}")
    tree_set = [n for n in keep if qc.loc[n, "source"] != "public" or n in selected]
    write_fasta(os.path.join(OUTDIR, "03_qc_cox1.fasta"), [(n, seqs[n]) for n in tree_set])
    ts = qc.loc[tree_set, ["source", "length", "ref_identity"]]
    ts["represents"] = [";".join(m for m in qc.index[qc.kept_as == n] if m != n)
                        for n in tree_set]
    ts.to_csv(os.path.join(OUTDIR, "03_qc_tree_set.tsv"), sep="\t")
    print(f"[03] tree set: {len(tree_set)} sequences "
          + str(ts.source.value_counts().to_dict()))


if __name__ == "__main__":
    sys.exit(main())
