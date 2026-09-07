"""Pruebas del núcleo científico (algoritmos puros)."""

import pytest

from app.core import alignment, fasta, restriction, sequences, tools


# ---------------------------------------------------------------- sequences
def test_reverse_complement_dna():
    assert sequences.reverse_complement("ATGC") == "GCAT"
    assert sequences.reverse_complement("GGATCC") == "GGATCC"  # palindrómico


def test_reverse_complement_rna():
    assert sequences.reverse_complement("AUGCAU", "rna") == "AUGCAU"


def test_gc_content():
    assert sequences.gc_content("GCGC") == 100.0
    assert sequences.gc_content("ATAT") == 0.0
    assert sequences.gc_content("AACC") == 50.0


def test_melting_temperature_reasonable_primers():
    tm = sequences.melting_temperature("ATGCTGACGTACGTCAGCTAGCTA")
    assert 40.0 < tm < 80.0


def test_composition_counts():
    comp = sequences.composition("AAGCT")
    assert comp == {"A": 2, "C": 1, "G": 1, "T": 1}


def test_invalid_symbols_raise():
    with pytest.raises(sequences.SequenceError):
        sequences.validate_nucleic("ATGX", "dna")  # X no es IUPAC


def test_empty_sequence_raises():
    with pytest.raises(sequences.SequenceError):
        sequences.validate_nucleic("   ")


# ------------------------------------------------------------------ tools
def test_translate_dna_to_protein():
    res = tools.translate("ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCC", frame=0)
    assert res["protein"].startswith("MVSKGEELFTGVVP")
    assert res["length_aa"] == len(res["protein"])


def test_translate_stop_codon():
    res = tools.translate("ATGGGCTAA", frame=0)
    assert res["protein"] == "MG*"
    assert res["has_stop_codon"] is True


def test_translate_frame_offsets():
    # Tres bases al inicio para desplazar el marco de lectura a 0.
    res = tools.translate("CGTATGGCGTAACGT", frame=0)
    assert res["protein"] == "RMA*"


def test_translate_accents_rna():
    res = tools.translate("AUGGGCUAG", kind="rna", frame=0)
    assert res["protein"] == "MG*"


def test_kmer_table():
    table = tools.kmer_table("AAAACCCC", 2)
    assert table["AA"] == 3
    assert table["CC"] == 3
    assert sum(table.values()) == 7


def test_find_orfs_both_strands():
    orfs = tools.find_orfs("ATGAAATAA", min_aa=1)
    assert orfs[0]["protein"] == "MK*"
    assert orfs[0]["strand"] == "+"


# ------------------------------------------------------------ restriction
def test_restriction_ecoRI_two_sites():
    cuts = restriction.restriction_sites("GAATTCAATTGAATTC")
    assert len(cuts["EcoRI"]["cuts"]) == 2


def test_restriction_palindrome():
    cuts = restriction.restriction_sites("GGATCC")
    c = cuts["BamHI"]["cuts"][0]
    # corte simétrico sobre la hebra complementaria
    assert c["sense_cut"] == 1


def test_restriction_unknown_enzyme():
    with pytest.raises(sequences.SequenceError):
        restriction.restriction_sites("ACGT", ["EnzimaFalsa"])


# --------------------------------------------------------------- alignment
def test_global_align_identity():
    result = alignment.align("GATTACAT", "GATTACAT")
    assert result["score"] == 16
    assert result["identity_percent"] == 100.0
    assert result["aligned_a"] == "GATTACAT"


def test_global_align_with_gap():
    result = alignment.align("GATTACA", "GATACA")
    # óptimo 10 (6 coincidencias × 2 − 1 gap × 2); el hueco puede ir
    # en cualquiera de las dos posiciones equivalentes.
    assert result["score"] == 10
    assert result["identity_percent"] == 85.71
    assert "-" in result["aligned_a"] or "-" in result["aligned_b"]
    assert len(result["aligned_a"]) == len(result["aligned_b"])


def test_local_align_finds_shared_motif():
    result = alignment.align("ACGTGGC", "GGCACGT", local=True)
    assert result["score"] > 4
    assert result["mode"] == "local"


def test_alignment_length_cap():
    with pytest.raises(ValueError):
        alignment.align("A" * 3000, "B" * 10)


# ------------------------------------------------------------------ fasta
def test_fasta_parse():
    records = fasta.parse_fasta(
        ">gfp enhanced green\nACGTACGT\n>seq2\nGGGGCCCC\n"
    )
    assert len(records) == 2
    assert records[0]["id"] == "gfp"
    assert records[0]["length"] == 8
    assert records[1]["gc_content"] == 100.0


def test_fasta_parse_multiline_and_spaces():
    records = fasta.parse_fasta(">a\nAC GT\nTT\n")
    assert records[0]["sequence"] == "ACGTTT"
    assert records[0]["length"] == 6


def test_fasta_empty_raises():
    with pytest.raises(sequences.SequenceError):
        fasta.parse_fasta("")