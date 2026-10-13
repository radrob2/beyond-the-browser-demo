#!/usr/bin/env python3
"""Step 01: download public cox1 for Hyalella + amphipod outgroups.

NCBI: direct E-utilities HTTP (esearch with history server -> efetch GenBank
flat), parsed with Biopython so source-feature metadata (country /
geo_loc_name, lat_lon, voucher, collection date) is kept.

BOLD (portal v5 API): all Hyalella COI-5P records. BOLD records mirrored in
GenBank (insdc_acs) are not duplicated; instead their BOLD metadata (country,
province, site, coordinates, BIN) fills gaps in the GenBank record. BOLD-only
records are added with their processid as header.

Outputs (in $WORK_DIR/01_download):
    01_download_hyalella_cox1_ncbi.fasta   GenBank Hyalella COI records
    01_download_hyalella_coi_bold.fasta    BOLD-only Hyalella COI-5P records
    01_download_outgroup_cox1_ncbi.fasta   one record per outgroup genus
    01_download_metadata.tsv               one row per downloaded record
"""
import io
import os
import sys
import time

import pandas as pd
import requests
from Bio import SeqIO

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API_KEY = os.environ.get("NCBI_API_KEY", "")
EMAIL = os.environ.get("NCBI_EMAIL", "")
OUTDIR = os.path.join(os.environ["WORK_DIR"], "01_download")

COI_TERMS = ('(COI[All Fields] OR COX1[All Fields] OR CO1[All Fields] '
             'OR "cytochrome oxidase subunit I"[All Fields] '
             'OR "cytochrome c oxidase subunit I"[All Fields] '
             'OR "cytochrome oxidase subunit 1"[All Fields] '
             'OR "cytochrome c oxidase subunit 1"[All Fields])')
HYALELLA_QUERY = f'Hyalella[Organism] AND {COI_TERMS} AND 300:2000[SLEN]'

# Outgroups: Hyalidae and Talitridae (fellow Talitroidea, closest to Hyalella),
# then increasingly distant amphipod families.
OUTGROUP_GENERA = [
    "Parhyale", "Hyale", "Allorchestes", "Talitrus", "Orchestia",
    "Platorchestia", "Gammarus", "Niphargus", "Caprella", "Ampithoe",
]
COX1_NAMES = {"cox1", "coi", "co1", "coxi", "cox-1", "coi-5p"}


def eutil(endpoint, **params):
    params.update({k: v for k, v in (("api_key", API_KEY), ("email", EMAIL)) if v})
    for attempt in range(5):
        r = requests.post(f"{EUTILS}/{endpoint}.fcgi", data=params, timeout=300)
        if r.status_code == 200:
            time.sleep(0.12 if API_KEY else 0.4)
            return r
        time.sleep(2 ** attempt)
    r.raise_for_status()


def search_and_fetch(term, retmax=100000):
    """esearch with usehistory, then efetch GenBank flat in batches."""
    r = eutil("esearch", db="nuccore", term=term, usehistory="y", retmax=0)
    txt = r.text
    count = int(txt.split("<Count>")[1].split("</Count>")[0])
    webenv = txt.split("<WebEnv>")[1].split("</WebEnv>")[0]
    qkey = txt.split("<QueryKey>")[1].split("</QueryKey>")[0]
    count = min(count, retmax)
    records = []
    for start in range(0, count, 500):
        r = eutil("efetch", db="nuccore", WebEnv=webenv, query_key=qkey,
                  rettype="gb", retmode="text", retstart=start,
                  retmax=min(500, count - start))
        records.extend(SeqIO.parse(io.StringIO(r.text), "genbank"))
    return records


def cox1_region(rec):
    """Return the cox1 sequence: the annotated CDS/gene if the record carries
    other loci too, otherwise the whole record."""
    for ftype in ("CDS", "gene"):
        for f in rec.features:
            if f.type != ftype:
                continue
            names = {v.lower() for k in ("gene", "product")
                     for v in f.qualifiers.get(k, [])}
            if names & COX1_NAMES or any("oxidase subunit i" in n or
                                         "oxidase subunit 1" in n for n in names):
                sub = f.extract(rec.seq)
                if len(sub) < len(rec.seq) - 50:
                    return str(sub)
                return str(rec.seq)
    return str(rec.seq)


def metadata(rec, group):
    src = next((f for f in rec.features if f.type == "source"), None)
    q = src.qualifiers if src else {}
    get = lambda k: q.get(k, [""])[0]
    country = get("geo_loc_name") or get("country")
    study = next((r.title for r in rec.annotations.get("references", [])
                  if r.title and r.title != "Direct Submission"), "")
    return {
        "accession": rec.id,
        "db": "ncbi",
        "group": group,
        "organism": rec.annotations.get("organism", get("organism")),
        "bin": "",
        "country": country.split(":")[0].strip(),
        "locality": country.split(":", 1)[1].strip() if ":" in country else "",
        "lat_lon": get("lat_lon"),
        "isolate": get("isolate") or get("specimen_voucher") or get("strain"),
        "collection_date": get("collection_date"),
        "study": study,
        "metadata_from_bold": False,
    }


def bold_fetch():
    """All Hyalella COI-5P records from the BOLD v5 portal API."""
    base = "https://portal.boldsystems.org/api"
    r = requests.get(f"{base}/query", params={"query": "tax:genus:Hyalella",
                                              "extent": "full"}, timeout=300)
    r.raise_for_status()
    qid = r.json()["query_id"]
    r = requests.get(f"{base}/documents/{qid}/download", params={"format": "tsv"},
                     timeout=600)
    r.raise_for_status()
    df = pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str,
                     keep_default_na=False, quoting=3)
    df = df[(df.marker_code == "COI-5P") & (df.nuc != "")].copy()
    df["nuc"] = df.nuc.str.replace("-", "", regex=False).str.upper()
    # BOLD placeholders -> empty; BOLD country names -> GenBank spelling
    df["country/ocean"] = df["country/ocean"].replace(
        {"Unrecoverable": "", "Unspecified country": "", "United States": "USA"})
    return df


def bold_locality(row):
    parts = [row["province/state"], row["region"], row["sector"], row["site"]]
    return ", ".join(dict.fromkeys(p for p in parts if p))


def write(records, group, fasta_path, rows):
    with open(fasta_path, "w") as fh:
        for rec in records:
            seq = cox1_region(rec).upper()
            fh.write(f">{rec.id} {rec.description}\n{seq}\n")
            row = metadata(rec, group)
            row["length"] = len(seq)
            rows.append(row)


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    rows = []

    print(f"[01] query: {HYALELLA_QUERY}")
    hy = search_and_fetch(HYALELLA_QUERY)
    print(f"[01] Hyalella records fetched: {len(hy)}")
    write(hy, "hyalella",
          os.path.join(OUTDIR, "01_download_hyalella_cox1_ncbi.fasta"), rows)

    outgroups = []
    for genus in OUTGROUP_GENERA:
        recs = search_and_fetch(f"{genus}[Organism] AND {COI_TERMS} AND 550:1600[SLEN]",
                                retmax=20)
        recs = [r for r in recs if "Hyalella" not in r.annotations.get("organism", "")]
        if not recs:
            print(f"[01] WARNING: no outgroup record for {genus}")
            continue
        # prefer a record with country metadata, then the longest
        recs.sort(key=lambda r: (not metadata(r, "")["country"], -len(r.seq)))
        outgroups.append(recs[0])
        print(f"[01] outgroup {genus}: {recs[0].id} "
              f"{recs[0].annotations.get('organism')}")
    write(outgroups, "outgroup",
          os.path.join(OUTDIR, "01_download_outgroup_cox1_ncbi.fasta"), rows)
    df = pd.DataFrame(rows)

    print("[01] BOLD: tax:genus:Hyalella, COI-5P")
    bold = bold_fetch()
    print(f"[01] BOLD COI-5P records with sequence: {len(bold)}")
    ncbi_base = {a.split(".")[0]: i for i, a in enumerate(df.accession)}
    bold["insdc_base"] = bold.insdc_acs.str.split(".").str[0]
    mirrored = bold[bold.insdc_base.isin(ncbi_base)]
    filled = 0
    for _, b in mirrored.iterrows():
        i = df.index[ncbi_base[b.insdc_base]]
        df.at[i, "bin"] = b.bin_uri
        if not df.at[i, "country"] and b["country/ocean"]:
            df.at[i, "country"] = b["country/ocean"]
            df.at[i, "locality"] = bold_locality(b)
            df.at[i, "lat_lon"] = df.at[i, "lat_lon"] or b.coord
            df.at[i, "metadata_from_bold"] = True
            filled += 1
        elif not df.at[i, "locality"] and b["country/ocean"] == df.at[i, "country"]:
            df.at[i, "locality"] = bold_locality(b)
    print(f"[01] BOLD records mirrored in GenBank: {len(mirrored)}; "
          f"GenBank records given country/locality from BOLD: {filled}")
    only = bold[~bold.insdc_base.isin(ncbi_base)].drop_duplicates("processid")
    with open(os.path.join(OUTDIR, "01_download_hyalella_coi_bold.fasta"), "w") as fh:
        for _, b in only.iterrows():
            fh.write(f">{b.processid}\n{b.nuc}\n")
    bold_rows = pd.DataFrame({
        "accession": only.processid, "db": "bold", "group": "hyalella",
        "organism": only.identification.where(only.identification != "", "Hyalella"),
        "bin": only.bin_uri,
        "country": only["country/ocean"],
        "locality": [bold_locality(r) for _, r in only.iterrows()],
        "lat_lon": only.coord, "isolate": only.sampleid,
        "collection_date": only.collection_date_start, "study": "",
        "metadata_from_bold": True, "length": only.nuc.str.len()})
    print(f"[01] BOLD-only records added: {len(bold_rows)}")
    df = pd.concat([df, bold_rows], ignore_index=True)

    df.to_csv(os.path.join(OUTDIR, "01_download_metadata.tsv"), sep="\t", index=False)
    hyd = df[df.group == "hyalella"]
    print(f"[01] Hyalella total: {len(hyd)} records ({(hyd.db == 'ncbi').sum()} GenBank, "
          f"{(hyd.db == 'bold').sum()} BOLD-only), {hyd.organism.nunique()} names, "
          f"{hyd.bin.replace('', None).nunique()} BINs, "
          f"{(hyd.country != '').sum()} with country")
    print(f"[01] outgroups: {(df.group == 'outgroup').sum()}")


if __name__ == "__main__":
    sys.exit(main())
