#!/usr/bin/env python3
"""Step 1: download public cox1 (COI) sequences + source metadata from NCBI.

Direct HTTP to NCBI E-utilities (stdlib urllib only; no entrez-direct).
Pattern: esearch (usehistory=y) -> efetch via WebEnv/query_key in batches.

Outputs (in --outdir):
    01_download_hyalella_cox1.fasta   every ingroup cox1 record in the length window
    01_download_outgroup_cox1.fasta   one cox1 sequence per outgroup amphipod genus
    01_download_public_metadata.tsv   one row per downloaded record (both sets)

Usage:
    01_download_public.py --outdir DIR --ingroup Hyalella --outgroup-n 40 \
        --slen-min 300 --slen-max 2000
"""
import argparse
import collections
import concurrent.futures
import json
import os
import re
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
API_KEY = os.environ.get("NCBI_API_KEY", "")
# NCBI allows 3 requests/sec without a key, 10/sec with one.
MIN_INTERVAL = 0.11 if API_KEY else 0.40
MAX_WORKERS = 6 if API_KEY else 3
MAX_TRIES = 7

FASTA_BATCH = 500
GB_BATCH = 200
OUTGROUP_MIN_LEN = 600      # preferred minimum; falls back to --slen-min
OUTGROUP_CANDIDATES_PER_GENUS = 40
RELATIVE_FRACTION = 0.6     # share of outgroup slots for Talitroidea/Hyaloidea

# cox1 in the gene field or spelled out in the title.
COX1_CLAUSE = (
    '(COI[Gene] OR COX1[Gene] OR CO1[Gene] OR COXI[Gene]'
    ' OR "cytochrome c oxidase subunit I"[Title]'
    ' OR "cytochrome oxidase subunit I"[Title]'
    ' OR "cytochrome oxidase subunit 1"[Title]'
    ' OR "cytochrome c oxidase subunit 1"[Title])'
)

# Ordered, hard-coded outgroup genera. Relatives of Hyalella first
# (Talitroidea / Hyaloidea: Hyalidae, Talitridae, Chiltoniidae,
# Dogielinotidae, Najnidae), then progressively more distant amphipods.
RELATIVE_GENERA = [
    "Hyale", "Parhyale", "Apohyale", "Ptilohyale", "Protohyale", "Serejohyale",
    "Platorchestia", "Orchestia", "Talitrus", "Talorchestia",
    "Chiltonia", "Austrochiltonia", "Allorchestes", "Parallorchestes",
    "Dogielinotus", "Transorchestia", "Cryptorchestia", "Bellorchestia",
    "Megalorchestia", "Americorchestia", "Sardorchestia", "Deshayesorchestia",
    "Britorchestia", "Speziorchestia", "Macarorchestia", "Arcitalitrus",
    "Talitroides", "Phreatochiltonia", "Najna", "Parhyalella", "Hyachelia",
    "Afrochiltonia", "Wangiannachiltonia", "Arabunnachiltonia",
]
DISTANT_GENERA = [
    "Gammarus", "Echinogammarus", "Niphargus", "Crangonyx", "Eulimnogammarus",
    "Dikerogammarus", "Pontogammarus", "Corophium", "Monoporeia", "Caprella",
    "Jassa", "Ampithoe", "Melita", "Bathyporeia", "Gammaracanthus",
    "Pallaseopsis", "Stygobromus", "Synurella", "Pontoporeia", "Leucothoe",
    "Eusirus", "Orchomenella",
]

SOUTH_AMERICA = ["Chile", "Peru", "Bolivia", "Argentina", "Brazil", "Uruguay",
                 "Colombia", "Ecuador", "Venezuela", "Paraguay"]

META_COLUMNS = ["accession", "set", "organism", "country", "locality",
                "lat_lon", "isolate", "voucher", "length", "title"]

# Titles that indicate something other than a plain cox1 fragment.
NOT_COX1_ONLY = re.compile(
    r"subunit (II|III|2|3)\b|\bCOII\b|\bCOIII\b|\bCOX2\b|\bCOX3\b|\btRNA\b"
    r"|ribosomal|complete genome|UNVERIFIED|pseudogene|\bnuclear\b", re.I)
UNNAMED_SPECIES = re.compile(r"\b(sp|cf|aff|nr|gen|complex)\b\.?|[0-9]", re.I)


def log(msg):
    print(msg, flush=True)


# ---------------------------------------------------------------- HTTP ----

_rate_lock = threading.Lock()
_next_slot = [0.0]
_n_requests = [0]


def _wait_for_slot():
    """Global rate limiter: space request starts MIN_INTERVAL apart."""
    with _rate_lock:
        now = time.monotonic()
        slot = max(now, _next_slot[0])
        _next_slot[0] = slot + MIN_INTERVAL
        _n_requests[0] += 1
    delay = slot - time.monotonic()
    if delay > 0:
        time.sleep(delay)


def eutils(endpoint, post=False, **params):
    """Call an E-utilities endpoint; retry 429/5xx/network errors with backoff."""
    if API_KEY:
        params["api_key"] = API_KEY
    query = urllib.parse.urlencode(params)
    url = EUTILS + endpoint
    for attempt in range(1, MAX_TRIES + 1):
        _wait_for_slot()
        try:
            if post:
                req = urllib.request.Request(url, data=query.encode())
            else:
                req = urllib.request.Request(url + "?" + query)
            with urllib.request.urlopen(req, timeout=120) as resp:
                return resp.read().decode("utf-8")
        except urllib.error.HTTPError as err:
            if err.code != 429 and err.code < 500:
                raise
            reason = "HTTP %d" % err.code
        except (urllib.error.URLError, TimeoutError, ConnectionError) as err:
            reason = repr(err)
        if attempt == MAX_TRIES:
            sys.exit("ERROR: %s failed after %d tries (%s)"
                     % (endpoint, MAX_TRIES, reason))
        backoff = min(2 ** attempt, 30)
        log("  retry %d/%d for %s after %s; sleeping %ds"
            % (attempt, MAX_TRIES - 1, endpoint, reason, backoff))
        time.sleep(backoff)


def esearch(term, retmax=0, usehistory=False):
    """Return (count, idlist, webenv, query_key) for a nuccore search."""
    params = dict(db="nuccore", term=term, retmode="json", retmax=retmax)
    if usehistory:
        params["usehistory"] = "y"
    for attempt in range(MAX_TRIES):
        res = json.loads(eutils("esearch.fcgi", **params))["esearchresult"]
        if "ERROR" not in res and "count" in res:
            break
        log("  esearch returned an error payload; retrying")
        time.sleep(2)
    else:
        sys.exit("ERROR: esearch kept failing for: " + term)
    return (int(res["count"]), res.get("idlist", []),
            res.get("webenv", ""), res.get("querykey", ""))


# ------------------------------------------------------------- parsing ----

def parse_fasta(text):
    """Return list of (header_without_gt, UPPERCASE single-line sequence)."""
    records = []
    header, parts = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(parts).upper()))
            header, parts = line[1:].rstrip(), []
        elif line.strip():
            parts.append(line.strip())
    if header is not None:
        records.append((header, "".join(parts).upper()))
    return records


def clean(value):
    """Collapse tabs/newlines/runs of whitespace; None -> ''."""
    return re.sub(r"\s+", " ", value or "").strip()


def parse_genbank_xml(text):
    """Parse efetch GenBank XML into {versioned_accession: metadata dict}.

    rettype=gb&retmode=xml returns GBSeq_* tags; rettype=gbc returns the
    same structure as INSDSeq_*. Both are handled.
    """
    out = {}
    root = ET.fromstring(text)
    for p in ("GB", "INSD"):
        for seq in root.iter(p + "Seq"):
            acc = clean(seq.findtext(p + "Seq_accession-version"))
            quals = {}
            for feat in seq.iter(p + "Feature"):
                if feat.findtext(p + "Feature_key") != "source":
                    continue
                for q in feat.iter(p + "Qualifier"):
                    name = q.findtext(p + "Qualifier_name")
                    quals.setdefault(
                        name, clean(q.findtext(p + "Qualifier_value")))
                break
            geo = quals.get("geo_loc_name") or quals.get("country") or ""
            country, _, locality = geo.partition(":")
            out[acc] = {
                "accession": acc,
                "organism": clean(seq.findtext(p + "Seq_organism")),
                "country": country.strip(),
                "locality": locality.strip(),
                "lat_lon": quals.get("lat_lon", ""),
                "isolate": quals.get("isolate", ""),
                "voucher": quals.get("specimen_voucher", ""),
                "length": clean(seq.findtext(p + "Seq_length")),
                "title": clean(seq.findtext(p + "Seq_definition")),
                "taxonomy": clean(seq.findtext(p + "Seq_taxonomy")),
            }
    return out


# ------------------------------------------------------------ fetching ----

def fetch_set(pool, label, expected, **source):
    """Fetch FASTA + GenBank XML for a record set, batches run concurrently.

    `source` is either WebEnv/query_key (history server) or id=<csv>.
    Returns (fasta_records, metadata_by_accession).
    """
    post = "id" in source
    jobs = []
    for start in range(0, expected, FASTA_BATCH):
        jobs.append(("fasta", pool.submit(
            eutils, "efetch.fcgi", post=post, db="nuccore", rettype="fasta",
            retmode="text", retstart=start, retmax=FASTA_BATCH, **source)))
    for start in range(0, expected, GB_BATCH):
        jobs.append(("gb", pool.submit(
            eutils, "efetch.fcgi", post=post, db="nuccore", rettype="gb",
            retmode="xml", retstart=start, retmax=GB_BATCH, **source)))
    log("  [%s] %d efetch batches queued (%d fasta, %d genbank-xml)"
        % (label, len(jobs), sum(k == "fasta" for k, _ in jobs),
           sum(k == "gb" for k, _ in jobs)))
    fasta, meta = [], {}
    for kind, job in jobs:
        text = job.result()
        if kind == "fasta":
            fasta.extend(parse_fasta(text))
        else:
            meta.update(parse_genbank_xml(text))
    log("  [%s] got %d fasta records, %d genbank records (expected %d)"
        % (label, len(fasta), len(meta), expected))
    if len(fasta) != expected or len(meta) != expected:
        sys.exit("ERROR: [%s] incomplete download; rerun" % label)
    return fasta, meta


def pick_outgroups(pool, n_wanted, slen_min, slen_max):
    """Choose one cox1 accession per genus. Returns list of (genus, group, uid, summary)."""
    genera = ([(g, "relative") for g in RELATIVE_GENERA]
              + [(g, "distant") for g in DISTANT_GENERA])

    def term_for(genus, lo):
        return ("%s[Organism] AND Amphipoda[Organism] AND %s AND %d:%d[SLEN]"
                % (genus, COX1_CLAUSE, lo, slen_max))

    def search(genus, lo):
        return esearch(term_for(genus, lo), retmax=OUTGROUP_CANDIDATES_PER_GENUS)

    lo_pref = max(OUTGROUP_MIN_LEN, slen_min)
    log("  outgroup query template: " + term_for("<Genus>", lo_pref))
    futures = {g: pool.submit(search, g, lo_pref) for g, _ in genera}
    hits = {g: f.result() for g, f in futures.items()}
    # Fallback to the shorter length window for genera with no long record.
    if lo_pref > slen_min:
        retry = [g for g, _ in genera if hits[g][0] == 0]
        futures = {g: pool.submit(search, g, slen_min) for g in retry}
        for g, f in futures.items():
            hits[g] = f.result()

    all_ids = sorted({i for h in hits.values() for i in h[1]}, key=int)
    summaries = {}
    jobs = [pool.submit(eutils, "esummary.fcgi", post=True, db="nuccore",
                        retmode="json", id=",".join(all_ids[i:i + 400]))
            for i in range(0, len(all_ids), 400)]
    for job in jobs:
        result = json.loads(job.result())["result"]
        for uid in result.get("uids", []):
            summaries[uid] = result[uid]

    def acc_key(acc):
        # natural sort: letters prefix, then number, then version
        m = re.match(r"([A-Z_]+)(\d+)\.(\d+)", acc)
        return (m.group(1), int(m.group(2)), int(m.group(3))) if m else (acc, 0, 0)

    chosen = {}
    for genus, group in genera:
        cands = []
        for uid in hits[genus][1]:
            s = summaries.get(uid)
            if not s:
                continue
            title, organism = s.get("title", ""), s.get("organism", "")
            if NOT_COX1_ONLY.search(title):
                continue
            if not organism.startswith(genus + " "):
                continue
            epithet = organism[len(genus) + 1:]
            unnamed = bool(UNNAMED_SPECIES.search(epithet))
            short = int(s.get("slen", 0)) < OUTGROUP_MIN_LEN
            cands.append((short, unnamed, acc_key(s["accessionversion"]), uid))
        if cands:
            cands.sort()
            chosen[genus] = (group, cands[0][3])
        log("    %-20s %-8s total_hits=%-6d usable_of_first_%d=%-3d %s"
            % (genus, group, hits[genus][0], OUTGROUP_CANDIDATES_PER_GENUS,
               len(cands),
               summaries[cands[0][3]]["accessionversion"] if cands else "SKIPPED"))

    rel = [g for g in RELATIVE_GENERA if g in chosen]
    dist = [g for g in DISTANT_GENERA if g in chosen]
    n_rel = min(len(rel), max(int(round(n_wanted * RELATIVE_FRACTION)),
                              n_wanted - len(dist)))
    n_dist = min(len(dist), n_wanted - n_rel)
    picked = rel[:n_rel] + dist[:n_dist]
    log("  outgroup genera: %d relatives (Talitroidea/Hyaloidea) + %d distant"
        " = %d (wanted %d)" % (n_rel, n_dist, len(picked), n_wanted))
    return [(g, chosen[g][0], chosen[g][1], summaries[chosen[g][1]])
            for g in picked]


# -------------------------------------------------------------- output ----

def write_fasta(path, records):
    with open(path, "w") as out:
        for header, seq in records:
            out.write(">%s\n%s\n" % (header, seq))


def table(counter, top, indent="    "):
    for key, n in counter.most_common(top):
        log("%s%6d  %s" % (indent, n, key if key else "(none)"))
    if len(counter) > top:
        log("%s        ... %d more distinct values" % (indent, len(counter) - top))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--ingroup", required=True)
    ap.add_argument("--outgroup-n", type=int, required=True)
    ap.add_argument("--slen-min", type=int, required=True)
    ap.add_argument("--slen-max", type=int, required=True)
    args = ap.parse_args()

    t0 = time.monotonic()
    hy_fasta_path = os.path.join(args.outdir, "01_download_hyalella_cox1.fasta")
    og_fasta_path = os.path.join(args.outdir, "01_download_outgroup_cox1.fasta")
    meta_path = os.path.join(args.outdir, "01_download_public_metadata.tsv")

    log("NCBI_API_KEY: %s -> %.2fs between requests, %d workers"
        % ("set" if API_KEY else "not set", MIN_INTERVAL, MAX_WORKERS))

    with concurrent.futures.ThreadPoolExecutor(MAX_WORKERS) as pool:
        # ---- ingroup -----------------------------------------------------
        hy_term = "%s[Organism] AND %s AND %d:%d[SLEN]" % (
            args.ingroup, COX1_CLAUSE, args.slen_min, args.slen_max)
        log("[1/3] ingroup esearch: " + hy_term)
        count, _, webenv, qkey = esearch(hy_term, usehistory=True)
        log("  %d records match" % count)
        if count == 0:
            sys.exit("ERROR: ingroup query returned nothing")
        hy_fasta, hy_meta = fetch_set(pool, "hyalella", count,
                                      WebEnv=webenv, query_key=qkey)

        # ---- outgroup ----------------------------------------------------
        log("[2/3] outgroup: one cox1 record per genus (%d genera tried)"
            % (len(RELATIVE_GENERA) + len(DISTANT_GENERA)))
        picks = pick_outgroups(pool, args.outgroup_n, args.slen_min,
                               args.slen_max)
        og_ids = ",".join(p[2] for p in picks)
        og_fasta, og_meta = fetch_set(pool, "outgroup", len(picks), id=og_ids)

    # ---- write (sorted by accession => deterministic) ----------------------
    log("[3/3] writing outputs to " + args.outdir)
    hy_fasta.sort(key=lambda r: r[0].split()[0])
    og_fasta.sort(key=lambda r: r[0].split()[0])
    write_fasta(hy_fasta_path, hy_fasta)
    write_fasta(og_fasta_path, og_fasta)

    rows = []
    for set_name, fasta, meta in (("hyalella", hy_fasta, hy_meta),
                                  ("outgroup", og_fasta, og_meta)):
        for header, seq in fasta:
            acc = header.split()[0]
            if acc not in meta:
                sys.exit("ERROR: no GenBank metadata for FASTA record " + acc)
            row = dict(meta[acc], set=set_name)
            if int(row["length"]) != len(seq):
                sys.exit("ERROR: length mismatch for " + acc)
            rows.append(row)
    with open(meta_path, "w") as out:
        out.write("\t".join(META_COLUMNS) + "\n")
        for row in rows:
            out.write("\t".join(str(row[c]) for c in META_COLUMNS) + "\n")

    # ---- sanity checks -----------------------------------------------------
    log("")
    log("=== sanity checks ===")
    hy_rows = [r for r in rows if r["set"] == "hyalella"]
    og_rows = [r for r in rows if r["set"] == "outgroup"]
    n_fasta = len(hy_fasta) + len(og_fasta)
    accs = [r["accession"] for r in rows]
    problems = []
    log("records: hyalella=%d  outgroup=%d  total=%d"
        % (len(hy_fasta), len(og_fasta), n_fasta))
    log("fasta records (%d) == metadata rows (%d): %s"
        % (n_fasta, len(rows), "OK" if n_fasta == len(rows) else "FAIL"))
    if n_fasta != len(rows):
        problems.append("fasta/metadata count mismatch")
    n_unique = len(set(accs))
    log("unique accessions: %d of %d: %s"
        % (n_unique, len(accs), "OK" if n_unique == len(accs) else "FAIL"))
    if n_unique != len(accs):
        problems.append("duplicate accessions")
    odd = [a for a in accs if not re.match(r"^[A-Z]{1,6}_?\d+\.\d+$", a)]
    log("headers not matching >ACCESSION.VERSION: %d %s" % (len(odd), odd[:5]))

    lens = sorted(int(r["length"]) for r in hy_rows)
    log("hyalella length: min=%d median=%d max=%d; <600bp: %d"
        % (lens[0], lens[len(lens) // 2], lens[-1],
           sum(1 for x in lens if x < 600)))
    not_ingroup = [r["accession"] for r in hy_rows
                   if not r["organism"].startswith(args.ingroup)]
    log("hyalella records whose organism is not %s*: %d %s"
        % (args.ingroup, len(not_ingroup), not_ingroup[:5]))
    suspicious = [r for r in hy_rows + og_rows
                  if NOT_COX1_ONLY.search(r["title"])]
    log("titles suggesting non-cox1 / multi-gene / unverified records: %d"
        % len(suspicious))
    for r in suspicious[:15]:
        log("    %s  %s" % (r["accession"], r["title"]))
    not_amphipod = [r["accession"] for r in og_rows
                    if "Amphipoda" not in meta_tax(og_meta, r)]
    log("outgroup records outside Amphipoda: %d %s"
        % (len(not_amphipod), not_amphipod))
    if not_amphipod:
        problems.append("non-amphipod outgroup")

    for field in ("country", "locality", "lat_lon", "isolate", "voucher"):
        log("hyalella records with %-8s: %4d of %d"
            % (field, sum(1 for r in hy_rows if r[field]), len(hy_rows)))
    countries = collections.Counter(r["country"] for r in hy_rows)
    log("hyalella records per country (top 25):")
    table(countries, 25)
    log("hyalella South American records: %d (%s)"
        % (sum(countries[c] for c in SOUTH_AMERICA),
           ", ".join("%s=%d" % (c, countries[c]) for c in SOUTH_AMERICA)))
    log("hyalella records per organism (top 30):")
    table(collections.Counter(r["organism"] for r in hy_rows), 30)

    log("outgroup records (genus-list order):")
    by_acc = {r["accession"]: r for r in og_rows}
    for genus, group, _, summ in picks:
        r = by_acc[summ["accessionversion"]]
        fam = family_of(og_meta[r["accession"]]["taxonomy"])
        log("    %-8s %-12s %5sbp  %-18s %-36s %s"
            % (group, r["accession"], r["length"], fam, r["organism"],
               r["country"]))
    og_short = sum(1 for r in og_rows if int(r["length"]) < OUTGROUP_MIN_LEN)
    log("outgroup records <%dbp: %d" % (OUTGROUP_MIN_LEN, og_short))

    log("")
    log("wrote %s (%d records)" % (hy_fasta_path, len(hy_fasta)))
    log("wrote %s (%d records)" % (og_fasta_path, len(og_fasta)))
    log("wrote %s (%d rows)" % (meta_path, len(rows)))
    log("E-utilities requests: %d; python wall time: %.1fs"
        % (_n_requests[0], time.monotonic() - t0))
    if problems:
        sys.exit("ERROR: sanity checks failed: " + "; ".join(problems))


def meta_tax(meta, row):
    return meta[row["accession"]]["taxonomy"]


def family_of(taxonomy):
    """Last lineage element ending in -idae (the family), else ''."""
    fams = [t.strip() for t in taxonomy.split(";") if t.strip().endswith("idae")]
    return fams[-1] if fams else ""


if __name__ == "__main__":
    main()
