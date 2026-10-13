#!/usr/bin/env python3
"""QC helper for 03_qc.sh.

  03_qc.py length IN OUTDIR MIN_LENGTH   -> 03_qc_length_cox1.fasta
  03_qc.py finish IN OUTDIR REFERENCE    -> 03_qc_cox1_predup.fasta (pre-dedup),
                                            03_qc_cox1.fasta (final),
                                            03_qc_status.tsv (one row per input)
"""
import collections
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cox1_common import classify, read_fasta, revcomp, seq_id, write_fasta

SOURCE_RANK = {'pcr': 0, 'genome': 1, 'public': 2}
SAFE_ID = re.compile(r'^[A-Za-z0-9_.\-]+$')


def unique_ids(records):
    """Trim headers to the id token; suffix only true collisions."""
    seen, out, renamed = collections.Counter(), [], 0
    for header, seq in records:
        token = seq_id(header)
        seen[token] += 1
        if seen[token] > 1:
            token = f'{token}_dup{seen[token]}'
            renamed += 1
        if not SAFE_ID.match(token):
            sys.exit(f'[03] ERROR: header not tree-safe: {token}')
        out.append((token, seq))
    return out, renamed


def length_filter(in_path, outdir, min_len):
    records, renamed = unique_ids(read_fasta(in_path))
    keep = [(h, s) for h, s in records if len(s.replace('-', '')) >= min_len]
    write_fasta(os.path.join(outdir, '03_qc_length_cox1.fasta'), keep)
    print(f'[03]   input {len(records)}, headers made unique {renamed}, '
          f'kept {len(keep)}, dropped short {len(records) - len(keep)}')


def finish(in_path, outdir, ref_path):
    records, _ = unique_ids(read_fasta(in_path))
    ref_ids = {seq_id(h) for h, _ in read_fasta(ref_path)}
    if ref_ids & {h for h, _ in records}:
        sys.exit('[03] ERROR: QC reference found in merged dataset')
    passed_len = {h for h, _ in read_fasta(os.path.join(outdir, '03_qc_length_cox1.fasta'))}
    hits = {}
    with open(os.path.join(outdir, '03_qc_reference_hits.tsv')) as f:
        for line in f:
            q, ident, strand, alnlen = line.rstrip('\n').split('\t')
            hits[q] = (float(ident), strand, int(alnlen))

    status, oriented = {}, []
    for h, s in records:
        row = {'id': h, 'source': classify(h), 'length': len(s), 'ref_identity': '',
               'strand': '', 'status': 'kept', 'duplicate_of': ''}
        status[h] = row
        if h not in passed_len:
            row['status'] = 'dropped_short'
            continue
        if h not in hits:
            row['status'] = 'dropped_low_identity'
            continue
        ident, strand, _ = hits[h]
        row['ref_identity'], row['strand'] = f'{ident:.1f}', strand
        oriented.append((h, revcomp(s) if strand == '-' else s.upper()))
    write_fasta(os.path.join(outdir, '03_qc_cox1_predup.fasta'), oriented)

    # Deduplicate exact-identical sequences; our own samples are kept as the
    # representative so they always stay visible in the tree.
    by_seq = collections.defaultdict(list)
    for h, s in oriented:
        by_seq[s].append(h)
    final = []
    for s, ids in by_seq.items():
        ids.sort(key=lambda h: (SOURCE_RANK[classify(h)], h))
        for dup in ids[1:]:
            status[dup]['status'] = 'dropped_duplicate'
            status[dup]['duplicate_of'] = ids[0]
        final.append((ids[0], s))
    order = {h: i for i, (h, _) in enumerate(records)}
    final.sort(key=lambda r: order[r[0]])
    assert len({h for h, _ in final}) == len(final)
    write_fasta(os.path.join(outdir, '03_qc_cox1.fasta'), final)

    cols = ['id', 'source', 'length', 'ref_identity', 'strand', 'status', 'duplicate_of']
    with open(os.path.join(outdir, '03_qc_status.tsv'), 'w') as out:
        out.write('\t'.join(cols) + '\n')
        for row in status.values():
            out.write('\t'.join(str(row[c]) for c in cols) + '\n')

    table = collections.Counter((r['source'], r['status']) for r in status.values())
    statuses = ['kept', 'dropped_short', 'dropped_low_identity', 'dropped_duplicate']
    print('[03]   source   ' + '  '.join(statuses))
    for src in ('public', 'pcr', 'genome'):
        print(f'[03]   {src:8s} ' + '  '.join(f'{table.get((src, st), 0):>{len(st)}d}' for st in statuses))
    flipped = sum(1 for r in status.values() if r['strand'] == '-')
    print(f'[03]   reverse-complemented to match reference strand: {flipped}')
    print(f'[03]   final unique sequences: {len(final)}')


if __name__ == '__main__':
    if sys.argv[1] == 'length':
        length_filter(sys.argv[2], sys.argv[3], int(sys.argv[4]))
    else:
        finish(sys.argv[2], sys.argv[3], sys.argv[4])
