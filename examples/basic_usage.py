#!/usr/bin/env python3
"""
examples/basic_usage.py
-----------------------
Demonstrates the mcr_ptm_annotator API.

Run from repo root:
    python examples/basic_usage.py
"""

from mcr_ptm_annotator import McrAPTMAnnotator, KNOWN_PTMS, report
from mcr_ptm_annotator.annotator import _MARBURGENSIS_MCRA_LENGTH

# ── 0. Print known PTM database ───────────────────────────────────────────
print("=" * 60)
print("Known MCR PTMs (from PDB 1MRO, Ermler et al. 1997)")
print("=" * 60)
for ptm in KNOWN_PTMS:
    print(f"  {ptm.name:<25} {ptm.residue_type}{ptm.position_marburgensis:<6} "
          f"[{ptm.pdb_code}]  {ptm.modification[:45]}")
print()

# ── 1. Annotate a synthetic sequence ─────────────────────────────────────
print("=" * 60)
print("Example 1: Synthetic sequence with PTM residues at correct positions")
print("=" * 60)

# Build a synthetic 553-aa sequence with correct residue types at PTM sites
seq = list("A" * _MARBURGENSIS_MCRA_LENGTH)
for ptm in KNOWN_PTMS:
    seq[ptm.position_marburgensis - 1] = ptm.residue_type
synthetic_seq = "".join(seq)

annotator = McrAPTMAnnotator(position_window=30)
hits = annotator.annotate_sequence(synthetic_seq, seq_id="synthetic_mcrA")
print(report.summary(hits, seq_id="synthetic_mcrA"))

# ── 2. Show PTM details ───────────────────────────────────────────────────
print("\n" + "=" * 60)
print("Example 2: Detailed hit information")
print("=" * 60)
for h in hits:
    print(f"\n  {h.ptm.name}")
    print(f"    Query position  : {h.query_position}")
    print(f"    Residue found   : {h.residue}")
    print(f"    Expected pos    : {h.expected_position}")
    print(f"    Delta           : {h.position_delta:+d}")
    print(f"    Confidence      : {h.confidence}")
    print(f"    Modification    : {h.ptm.modification}")
    print(f"    References      : {h.ptm.references[0]}")

# ── 3. Export ─────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("Example 3: Export to TSV and JSON")
print("=" * 60)

import json
tsv = report.to_tsv(hits, seq_id="synthetic_mcrA")
print("TSV (first 2 lines):")
for line in tsv.split("\n")[:2]:
    print(" ", line)

js = report.to_json(hits, seq_id="synthetic_mcrA")
parsed = json.loads(js)
print(f"\nJSON: {len(parsed)} hits serialised")

# ── 4. Scan real sequence (user must provide) ─────────────────────────────
print("\n" + "=" * 60)
print("Example 4: Scan a real McrA sequence from FASTA")
print("=" * 60)
print("NOTE: example_mcrA.fasta contains a synthetic placeholder.")
print("Replace with a real McrA sequence downloaded from:")
print("  UniProt Q8TN35 (M. acetivorans) or UniProt P11558 (M. marburgensis)")

results = annotator.annotate_fasta("examples/example_mcrA.fasta")
for seq_id, seq_hits in results.items():
    print(f"\n--- {seq_id} ---")
    print(report.summary(seq_hits, seq_id=seq_id))

print("\nDone.")
print("\n⚠ All positions are heuristic estimates. Confirm by alignment")
print("  against PDB 1MRO chain A and by mass spectrometry.")
