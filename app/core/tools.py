"""Traducción, detección de ORFs y análisis de k-mers.

Herramientas de biología molecular computacional:
  * traducción ADN/ARN → proteína (código genético estándar),
  * búsqueda de marcos de lectura abiertos (ORFs) en ambas hebras,
  * frecuencia de palabras cortas (k-mers) usada p. ej. en
    identificación de especies (huella genómica).
"""

from .sequences import (
    IUPAC_PROTEIN,
    SequenceError,
    normalize,
    reverse_complement,
    validate_nucleic,
)

# Código genético estándar (codones de ARN → aminoácido).
CODON_TABLE: dict[str, str] = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

DNA_TO_RNA = str.maketrans("T", "U")
STOP_CODONS = {"TAA", "TAG", "TGA"}


def _to_rna(sequence: str, kind: str) -> str:
    if kind == "rna":
        return sequence
    return sequence.translate(DNA_TO_RNA)


def translate(sequence: str, kind: str = "dna", frame: int = 0) -> dict:
    """Traduce un ORF a proteína.

    `frame` es el desplazamiento en lectura (0, 1 o 2) desde el extremo
    5' de la hebra sensora. Reporta la proteína y el codón de parada.
    """
    if frame not in (0, 1, 2):
        raise SequenceError("El frame debe ser 0, 1 o 2.")
    seq = validate_nucleic(sequence, kind)
    rna = _to_rna(seq, kind)
    codons = [rna[i : i + 3] for i in range(frame, len(rna) - 2, 3)]

    protein = []
    for codon in codons:
        aa = CODON_TABLE.get(codon)
        if aa is None:
            raise SequenceError(f"Codón no reconocido: {codon}")
        protein.append(aa)
        if aa == "*":
            break

    peptide = "".join(protein)
    return {
        "frame": frame,
        "length_aa": len(peptide),
        "protein": peptide,
        "has_stop_codon": peptide.endswith("*"),
    }


def find_orfs(sequence: str, kind: str = "dna", min_aa: int = 20) -> list[dict]:
    """Busca marcos de lectura abiertos en ambas hebras.

    Un ORF va desde el codón de inicio ATG hasta el primer codón de
    parada en el mismo marco. Es el paso previo a identificar una
    posible proteína en una secuencia desconocida.
    """
    if min_aa < 1:
        raise SequenceError("min_aa debe ser >= 1.")
    seq = validate_nucleic(sequence, kind)
    if kind == "rna":
        seq = seq.translate(str.maketrans("U", "T"))

    found: list[dict] = []
    for strand, target in (("+", seq), ("-", reverse_complement(seq))):
        for frame in range(3):
            # Escaneo en tríos (en fase): tras un ATG solo cuentan los
            # codones siguientes, nunca un desplazamiento de 1 base.
            start = None
            for i in range(frame, len(target) - 2, 3):
                codon = target[i : i + 3]
                if start is None:
                    if codon == "ATG":
                        start = i
                elif codon in STOP_CODONS:
                    orf_len = (i + 3 - start) // 3
                    if orf_len >= min_aa:
                        protein = "".join(
                            CODON_TABLE[_to_rna(target[j : j + 3], "dna")]
                            for j in range(start, i + 3, 3)
                        )
                        found.append({
                            "strand": strand,
                            "frame": frame,
                            "start": start,          # 0-based, hebra sensora
                            "end": i + 2,            # última base del stop
                            "length_bp": i + 3 - start,
                            "length_aa": orf_len,
                            "sequence": target[start : i + 3],
                            "protein": protein,
                        })
                    start = None

    found.sort(key=lambda o: o["length_bp"], reverse=True)
    return found


def canonical_protein(protein: str) -> str:
    """Valida una secuencia proteica IUPAC y la normaliza."""
    prot = normalize(protein)
    invalid = {c for c in prot if c not in IUPAC_PROTEIN}
    if invalid:
        raise SequenceError(f"Aminoácidos inválidos: {''.join(sorted(invalid))}")
    return prot


def kmer_table(sequence: str, k: int, kind: str = "dna") -> dict:
    """Frecuencia de k-mers (palabras de `k` bases) en la secuencia."""
    if not 1 <= k <= min(12, len(sequence)):
        raise SequenceError(f"k debe estar entre 1 y min(12, len(secuencia)).")
    seq = validate_nucleic(sequence, kind)
    counts: dict[str, int] = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        counts[kmer] = counts.get(kmer, 0) + 1
    return dict(sorted(counts.items(), key=lambda item: (-item[1], item[0])))