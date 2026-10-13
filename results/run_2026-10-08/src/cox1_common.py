#!/usr/bin/env python3
"""Shared helpers: FASTA I/O and source classification from header pattern."""
import re

ACCESSION_RE = re.compile(r'^[A-Z]{1,6}_?\d+\.\d+$')
COMP = str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN')


def read_fasta(path):
    """Return list of (header_without_gt, sequence)."""
    records, header, parts = [], None, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(parts)))
                header, parts = line[1:].rstrip('\r\n'), []
            else:
                parts.append(line.strip())
    if header is not None:
        records.append((header, ''.join(parts)))
    return records


def write_fasta(path, records):
    """Single-line, uppercase records (project FASTA convention)."""
    with open(path, 'w') as out:
        for header, seq in records:
            out.write(f'>{header}\n{seq.upper()}\n')


def seq_id(header):
    """Leading identifier token of a header (what tree builders keep)."""
    return header.split()[0]


def classify(header):
    """Source from header pattern alone: pcr / genome / public / unknown."""
    token = seq_id(header)
    if token.startswith('genome_'):
        return 'genome'
    if '_pcr_' in token or token.endswith('_pcr'):
        return 'pcr'
    if ACCESSION_RE.match(token):
        return 'public'
    return 'unknown'


def revcomp(seq):
    return seq.upper().translate(COMP)[::-1]


def read_tsv(path):
    with open(path) as f:
        cols = f.readline().rstrip('\n').split('\t')
        return [dict(zip(cols, line.rstrip('\n').split('\t'))) for line in f if line.strip()]
