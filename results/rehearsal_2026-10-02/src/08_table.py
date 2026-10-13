#!/usr/bin/env python3
"""Step 08: nearest public relative of each of our samples.

Uses the exhaustive vsearch search from step 03 (each of our sequences vs
ALL QC-passed public Hyalella, not only the tree subset). Identity =
matching bases / aligned columns, terminal gaps ignored.

Outputs:
    $SHARE_DIR/08_nearest_relatives.tsv          one row per sample of ours
    $SHARE_DIR/08_nearest_relatives_grouped.tsv  samples grouped by lineage
"""
import os
import sys

import pandas as pd

WORK = os.environ["WORK_DIR"]
SHARE = os.environ["SHARE_DIR"]
NN = os.path.join(WORK, "03_qc", "03_qc_nearest_public.tsv")
QC = os.path.join(WORK, "03_qc", "03_qc_report.tsv")
META = os.path.join(WORK, "01_download", "01_download_metadata.tsv")

# % identity bands (cox1 barcode; BOLD BINs are built around ~2 % divergence)
SAME = 98.0
BORDER = 96.0


def assess(pid):
    if pid >= SAME:
        return "matches a known lineage"
    if pid >= BORDER:
        return "close to a known lineage (borderline)"
    return "divergent: candidate new lineage"


# Context from the publication title when the record itself lacks it.
# Reported in its own, explicitly labelled column; never written into country.
STUDY_HINTS = [
    ("titicaca", "Lake Titicaca, Peru/Bolivia (from publication title)"),
    # Major et al.: lab strain vs wild strains from eastern North America
    ("laboratory strain", "eastern North America, lab or wild (from publication title)"),
    ("toxicity", "laboratory culture used in toxicity tests (from publication title)"),
]


def study_hint(title):
    t = title.lower()
    return next((hint for key, hint in STUDY_HINTS if key in t), "")


def short_location(country, locality):
    """'Canada, Ontario' style: country plus first locality element."""
    if not country:
        return ""
    first = locality.split(",")[0].strip() if locality else ""
    return f"{country}, {first}" if first else country


def species_label(r):
    org = r.organism or "Hyalella"
    if org == "Hyalella":
        org = "Hyalella sp. (unnamed)"
    return org


def main():
    meta = pd.read_csv(META, sep="\t", keep_default_na=False).set_index("accession")
    qc = pd.read_csv(QC, sep="\t", keep_default_na=False).set_index("name")
    nn = pd.read_csv(NN, sep="\t", keep_default_na=False)
    nn = nn[nn.target != "*"].sort_values(["query", "identity", "alnlen"],
                                          ascending=[True, False, False])
    ours = qc.index[qc.source.isin(["pcr", "genome"])]

    # names attached to each BIN anywhere in the public data
    named = meta[(meta.bin != "") & ~meta.organism.isin(["Hyalella", ""])]
    bin_names = named.groupby("bin").organism.agg(lambda s: "; ".join(sorted(set(s))))

    rows = []
    for q in ours:
        g = nn[nn["query"] == q]
        if g.empty:
            rows.append({"our_sample": q, "source": qc.loc[q, "source"],
                         "assessment": "no public hit >= 70 %"})
            continue
        top = g.iloc[0]
        m = meta.loc[top.target]
        tied = g[g.identity >= top.identity - 1e-9]
        # closest relative that does have a country, if the top hit lacks one
        geo = g[[bool(meta.loc[t, "country"]) for t in g.target]]
        geo_row = geo.iloc[0] if len(geo) else None
        chile = g[[meta.loc[t, "country"] == "Chile" for t in g.target]]
        rows.append({
            "our_sample": q,
            "source": qc.loc[q, "source"],
            "nearest_public_relative": top.target,
            "database": m.db,
            "relative_species": species_label(m),
            "relative_bin": m.bin,
            "bin_also_named_as": bin_names.get(m.bin, "") if m.bin else "",
            "relative_country": m.country or "not recorded",
            "relative_locality": m.locality or "not recorded",
            "pct_identity": round(top.identity, 1),
            "aligned_bp": int(top.alnlen),
            "n_equally_close": len(tied),
            "countries_of_equally_close": "; ".join(sorted(
                {meta.loc[t, "country"] for t in tied.target} - {""})),
            "nearest_with_country": "" if geo_row is None else geo_row.target,
            "nearest_with_country_species": "" if geo_row is None
                else species_label(meta.loc[geo_row.target]),
            "nearest_with_country_location": "" if geo_row is None else ", ".join(
                x for x in (meta.loc[geo_row.target, "country"],
                            meta.loc[geo_row.target, "locality"]) if x),
            "nearest_with_country_pct_identity": "" if geo_row is None
                else round(geo_row.identity, 1),
            # local context for an Atacama study: closest Chilean record, if
            # any is among the 25 closest public hits kept in step 03
            "closest_chilean_record": chile.iloc[0].target if len(chile)
                else "none among 25 closest",
            "closest_chilean_location": meta.loc[chile.iloc[0].target, "locality"]
                if len(chile) else "",
            "closest_chilean_pct_identity": round(chile.iloc[0].identity, 1)
                if len(chile) else "",
            "study_of_relative": m.study,
            "relative_context_from_study": study_hint(m.study),
            "assessment": assess(top.identity),
        })
    df = pd.DataFrame(rows).sort_values(["source", "our_sample"],
                                        ascending=[False, True])
    os.makedirs(SHARE, exist_ok=True)
    out = os.path.join(SHARE, "08_nearest_relatives.tsv")
    df.to_csv(out, sep="\t", index=False)
    print(f"[08] {len(df)} samples -> {out}")
    print(df.assessment.value_counts().to_string())

    # grouped view: one row per species of nearest relative
    def where(r):
        if r.relative_country != "not recorded":
            loc = short_location(r.relative_country, meta.loc[r.nearest_public_relative,
                                                              "locality"])
        elif r.relative_context_from_study:
            loc = r.relative_context_from_study.replace(" (from publication title)", "*")
        elif r.nearest_with_country_location:
            loc = short_location(*(r.nearest_with_country_location.split(", ", 1)
                                   + [""])[:2]) + "†"
        else:
            loc = "not recorded"
        return loc
    df["where"] = [where(r) for r in df.itertuples()]
    short = lambda n: n.split("_pcr_")[0] if "_pcr_" in n else \
        n.replace("genome_hyalella_sp_", "").replace("_cox1_extracted", "")
    grp = df.groupby("relative_species").agg(
        n_samples=("our_sample", "size"),
        n_pcr=("source", lambda s: (s == "pcr").sum()),
        n_genome=("source", lambda s: (s == "genome").sum()),
        samples=("our_sample", lambda s: ", ".join(short(x) for x in s)),
        min_pct_identity=("pct_identity", "min"),
        max_pct_identity=("pct_identity", "max"),
        relative_bins=("relative_bin", lambda s: ", ".join(sorted(set(s) - {""}))),
        closest_relative_locations=("where", lambda s: "; ".join(sorted(set(s)))),
        n_with_chilean_hit=("closest_chilean_pct_identity", lambda s: (s != "").sum()),
        chilean_min_pct=("closest_chilean_pct_identity",
                         lambda s: min([v for v in s if v != ""], default="")),
        chilean_max_pct=("closest_chilean_pct_identity",
                         lambda s: max([v for v in s if v != ""], default="")),
        chilean_regions=("closest_chilean_location", lambda s: "; ".join(
            sorted({v.split(",")[0] for v in s if v}))),
    ).sort_values("n_samples", ascending=False)
    # * = from publication title; † = nearest relative WITH a country (lower identity)
    gout = os.path.join(SHARE, "08_nearest_relatives_grouped.tsv")
    grp.to_csv(gout, sep="\t")
    print(f"[08] {len(grp)} lineage groups -> {gout}")


if __name__ == "__main__":
    sys.exit(main())
