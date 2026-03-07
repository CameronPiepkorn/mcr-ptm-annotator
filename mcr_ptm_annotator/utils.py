"""
mcr_ptm_annotator.utils
~~~~~~~~~~~~~~~~~~~~~~~
Shared utilities.
"""

from __future__ import annotations
from pathlib import Path
from typing import Iterator

# Standard single-letter amino acid codes
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def parse_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """
    Minimal FASTA parser for protein sequences.
    Yields (seq_id, sequence) tuples.
    """
    path = Path(path)
    current_id: str | None = None
    buffer: list[str] = []

    with path.open() as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                if current_id is not None:
                    yield current_id, "".join(buffer)
                current_id = line[1:].split()[0]
                buffer = []
            else:
                buffer.append(line.upper())

    if current_id is not None:
        yield current_id, "".join(buffer)


def is_protein_sequence(seq: str) -> bool:
    """Return True if sequence contains only standard amino acid characters."""
    return all(c in AMINO_ACIDS for c in seq.upper() if c not in (" ", "\n"))
