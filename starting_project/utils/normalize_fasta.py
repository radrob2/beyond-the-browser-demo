#!/usr/bin/env python3
"""Standardize a FASTA file: single-line sequences, uppercase letters.

Usage:
    ./utils/normalize_fasta.py INPUT [OUTPUT]

If OUTPUT is omitted, INPUT is rewritten in place via a tempfile.
"""
import os
import sys
import tempfile


def normalize(in_path, out_path):
    """Read FASTA, uppercase sequences, write single-line records."""
    records = []
    header = None
    seq_parts = []
    with open(in_path) as f:
        for line in f:
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(seq_parts)))
                header = line.rstrip('\r\n')
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None:
            records.append((header, ''.join(seq_parts)))

    with open(out_path, 'w') as out:
        for header, seq in records:
            out.write(header + '\n')
            out.write(seq.upper() + '\n')


def main(argv):
    if len(argv) < 2 or argv[1] in ('-h', '--help'):
        print(__doc__, file=sys.stderr)
        sys.exit(2 if len(argv) < 2 else 0)

    in_path = argv[1]
    out_path = argv[2] if len(argv) > 2 else None

    if out_path is None:
        in_dir = os.path.dirname(os.path.abspath(in_path)) or '.'
        with tempfile.NamedTemporaryFile('w', dir=in_dir, suffix='.fasta',
                                         delete=False) as tmp:
            tmp_path = tmp.name
        try:
            normalize(in_path, tmp_path)
            os.replace(tmp_path, in_path)
        except Exception:
            os.unlink(tmp_path)
            raise
    else:
        normalize(in_path, out_path)


if __name__ == '__main__':
    main(sys.argv)
