"""Utilidades sobre secuencias de ácidos nucleicos (ADN/ARN).

Cubre las operaciones básicas de un flujo de trabajo bioinformático:
validación IUPAC, composición de bases, contenido GC y temperatura
de melting (Tm) de primers/cDNA.
"""

IUPAC_DNA = set("ACGTRYSWKMBDHVN")
IUPAC_RNA = set("ACGURYSWKMBDHVN")
IUPAC_PROTEIN = set("ARNDCQEGHILKMFPSTWYVBZX*")

COMPLEMENT_DNA = {
    "A": "T", "T": "A", "G": "C", "C": "G",
    "R": "Y", "Y": "R", "S": "S", "W": "W",
    "K": "M", "M": "K", "B": "V", "V": "B",
    "D": "H", "H": "D", "N": "N",
}
COMPLEMENT_RNA = {"A": "U", "U": "A", "G": "C", "C": "G", "N": "N"}


class SequenceError(ValueError):
    """La secuencia o algún parámetro es inválido."""


def normalize(sequence: str) -> str:
    """Limpia espacios, normaliza a mayúsculas y valida que no esté vacía."""
    seq = "".join(sequence.split()).upper()
    if not seq:
        raise SequenceError("La secuencia está vacía.")
    return seq


def validate_nucleic(sequence: str, kind: str = "dna") -> str:
    """Valida una secuencia contra los códigos IUPAC y la normaliza."""
    seq = normalize(sequence)
    allowed = IUPAC_DNA if kind == "dna" else IUPAC_RNA
    invalid = {c for c in seq if c not in allowed}
    if invalid:
        raise SequenceError(
            f"Símbolos IUPAC inválidos para {kind.upper()}: "
            f"{''.join(sorted(invalid))}. "
            f"Permitidos: {''.join(sorted(allowed))}"
        )
    return seq


def composition(sequence: str, kind: str = "dna") -> dict[str, int]:
    """Conteo de cada base de la secuencia."""
    seq = validate_nucleic(sequence, kind)
    return {base: seq.count(base) for base in sorted(set(seq))}


def gc_content(sequence: str, kind: str = "dna") -> float:
    """Porcentaje de bases G+C (parámetro clave de estabilidad del ADN)."""
    seq = validate_nucleic(sequence, kind)
    gc = seq.count("G") + seq.count("C")
    return round(gc / len(seq) * 100, 2)


def reverse_complement(sequence: str, kind: str = "dna") -> str:
    """Cadena complementaria reversa (la que usa la polimerasa in vitro)."""
    seq = validate_nucleic(sequence, kind)
    table = COMPLEMENT_DNA if kind == "dna" else COMPLEMENT_RNA
    return "".join(table[c] for c in reversed(seq))


def melting_temperature(sequence: str, kind: str = "dna") -> float | None:
    """Tm por la fórmula de Marmur-Doty (ajustada a sal 50 mM).

    Estimación de la temperatura a la que se separan las hebras del ADN;
    esencial para diseñar primers de PCR. Devuelve None si la secuencia
    es demasiado corta (< 10 bases) para ser fiable.
    """
    seq = validate_nucleic(sequence, kind)
    gc = seq.count("G") + seq.count("C")
    n = len(seq)
    if n < 10:
        return None
    tm = 64.9 + 41 * (gc - 16.4) / n
    return round(max(0.0, tm), 2)