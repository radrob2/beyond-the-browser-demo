#!/usr/bin/env python3
"""Helper for 04_nearest_subset.sh: nearest-relative table and tree subset."""
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cox1_common import classify, read_fasta, read_tsv, write_fasta

work, share = sys.argv[2], sys.argv[3]
qc, out = os.path.join(work, '03_qc'), os.path.join(work, '04_subset')


def split():
    predup = read_fasta(os.path.join(qc, '03_qc_cox1_predup.fasta'))
    final = read_fasta(os.path.join(qc, '03_qc_cox1.fasta'))
    own = [(h, s) for h, s in predup if classify(h) != 'public']
    public = [(h, s) for h, s in final if classify(h) == 'public']
    write_fasta(os.path.join(out, '04_own_all.fasta'), own)
    write_fasta(os.path.join(out, '04_public_all.fasta'), public)
    print(f'[04] own samples passing QC: {len(own)}; public unique sequences: {len(public)}')


def finish():
    same, close, min_overlap = float(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6])
    k_near, per_species, n_out = int(sys.argv[7]), int(sys.argv[8]), int(sys.argv[9])
    meta = {r['accession']: r for r in read_tsv(os.path.join(work, '01_download', '01_download_public_metadata.tsv'))}
    status = {r['id']: r for r in read_tsv(os.path.join(qc, '03_qc_status.tsv'))}

    hits = collections.defaultdict(list)
    with open(os.path.join(out, '04_nearest_hits.tsv')) as f:
        for line in f:
            q, t, ident, alnlen, mism, gaps = line.rstrip('\n').split('\t')
            if int(alnlen) >= min_overlap:
                hits[q].append((float(ident), int(alnlen), t))
    for q in hits:
        hits[q].sort(key=lambda x: (-x[0], -x[1], x[2]))

    def verdict(ident):
        if ident >= same:
            return 'matches a known lineage'
        if ident >= close:
            return 'distinct, with a known close relative'
        return 'no close public match (candidate new lineage)'

    cols = ['our_sample', 'source', 'nearest_public_relative', 'relative_species',
            'relative_country', 'relative_locality', 'relative_lat_lon',
            'percent_identity', 'aligned_bp', 'n_public_equally_close',
            'countries_of_equally_close', 'interpretation', 'qc_status']
    rows, neighbours = [], []
    for sid, st in status.items():
        if st['source'] == 'public':
            continue
        row = dict.fromkeys(cols, '')
        row.update(our_sample=sid, source=st['source'], qc_status=st['status'])
        if st['status'] in ('kept', 'dropped_duplicate') and hits.get(sid):
            if st['status'] == 'dropped_duplicate':
                row['qc_status'] = f"identical to {st['duplicate_of']}"
            ident = hits[sid][0][0]
            tied = [t for i, _, t in hits[sid] if i == ident]
            # among equally close relatives, report the best-documented one
            alnlen, acc = min(((a, t) for i, a, t in hits[sid] if i == ident),
                              key=lambda x: (meta[x[1]]['country'] == '',
                                             meta[x[1]]['locality'] == '', -x[0], x[1]))
            m = meta[acc]
            tied_countries = sorted({meta[t]['country'] for t in tied if meta[t]['country']})
            row.update(nearest_public_relative=acc, relative_species=m['organism'],
                       relative_country=m['country'] or 'not recorded',
                       relative_locality=m['locality'] or 'not recorded',
                       relative_lat_lon=m['lat_lon'], percent_identity=f'{ident:.1f}',
                       aligned_bp=alnlen, n_public_equally_close=len(tied),
                       countries_of_equally_close='; '.join(tied_countries),
                       interpretation=verdict(ident))
            neighbours += [acc] + [t for _, _, t in hits[sid][:k_near]]
        else:
            row['interpretation'] = 'failed QC - not placed'
        rows.append(row)
    rows.sort(key=lambda r: (r['source'], r['our_sample']))
    with open(os.path.join(share, 'table', '04_nearest_public_relatives.tsv'), 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join(str(r[c]) for c in cols) + '\n')
    print(f'[04] table rows (one per sample of ours): {len(rows)}')
    for k, v in collections.Counter(r['interpretation'] for r in rows).most_common():
        print(f'[04]   {v:3d}  {k}')

    # ---- tree subset -----------------------------------------------------
    final = read_fasta(os.path.join(qc, '03_qc_cox1.fasta'))
    seqs = dict(final)
    keep = {h: 'own' for h in seqs if classify(h) != 'public'}
    for acc in neighbours:
        keep.setdefault(acc, 'neighbour')

    def rank(acc):  # prefer records with a country, then longer, then accession
        return (meta[acc]['country'] == '', -len(seqs[acc]), acc)

    by_species = collections.defaultdict(list)
    for acc in seqs:
        if acc in meta and meta[acc]['set'] == 'hyalella':
            by_species[meta[acc]['organism']].append(acc)
    for species, accs in sorted(by_species.items()):
        # backbone = formally named species only; provisional names
        # ("Hyalella sp. A1", "cf.", "aff.") enter only as nearest relatives
        if any(t in species.split() for t in ('sp.', 'cf.', 'aff.')):
            continue
        for acc in sorted(accs, key=rank)[:per_species]:
            keep.setdefault(acc, 'backbone')

    # outgroups: the n_out most reference-like genera (shortest branches)
    by_genus = {}
    outs = [a for a in seqs if a in meta and meta[a]['set'] == 'outgroup']
    for acc in sorted(outs, key=lambda a: (-float(status[a]['ref_identity']), a)):
        by_genus.setdefault(meta[acc]['organism'].split()[0], acc)
    for acc in list(by_genus.values())[:n_out]:
        keep[acc] = 'outgroup'

    subset = [(h, s) for h, s in final if h in keep]
    write_fasta(os.path.join(out, '04_subset_cox1.fasta'), subset)
    with open(os.path.join(out, '04_subset_roles.tsv'), 'w') as f:
        f.write('id\trole\n')
        for h, _ in subset:
            f.write(f'{h}\t{keep[h]}\n')
    roles = collections.Counter(keep.values())
    print(f'[04] tree subset: {len(subset)} of {len(final)} QC-passed sequences '
          f'({", ".join(f"{k} {v}" for k, v in sorted(roles.items()))})')
    print(f'[04] outgroups passing QC: {len(outs)} sequences / {len(by_genus)} genera; kept {roles.get("outgroup", 0)}')


split() if sys.argv[1] == 'split' else finish()
