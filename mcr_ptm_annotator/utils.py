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


DNA_CHARS = set("ATCGNU")  # nucleotide characters (RNA U included)


def is_protein_sequence(seq: str) -> bool:
    """
    Return True if sequence looks like a protein (not DNA/RNA).

    A sequence is treated as DNA/RNA — and therefore NOT a protein — if
    all of its characters fall within the nucleotide alphabet (ATCGNU).
    This handles the ambiguity that A, C, G, T are also valid amino acid
    single-letter codes.
    """
    cleaned = seq.upper().replace(" ", "").replace("\n", "")
    if not cleaned:
        return False
    # Reject if every character is a nucleotide symbol
    if all(c in DNA_CHARS for c in cleaned):
        return False
    return all(c in AMINO_ACIDS for c in cleaned)
