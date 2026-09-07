"""Parseo del formato FASTA (formato estándar de secuencias biológicas).

El FASTA es el formato de texto que usan GenBank/NCBI y los institutos
de investigación para compartir secuencias de ADN/proteínas. Cada
registro empieza con ">" seguido del identificador y una descripción.
"""

from .sequences import SequenceError, gc_content, validate_nucleic


def _kind_of(seq: str) -> str:
    """Infiera ADN/ARN: presencia de U sin T sugiere ARN."""
    return "rna" if "U" in seq and "T" not in seq else "dna"


def parse_fasta(content: str) -> list[dict]:
    """Convierte texto FASTA en una lista de registros analizados."""
    records: list[dict] = []
    current_id = current_desc = ""
    current_seq: list[str] = []

    def flush() -> None:
        seq = "".join(current_seq).upper().replace(" ", "")
        if seq:
            kind = _kind_of(seq)
            validate_nucleic(seq, kind)
            records.append({
                "id": (current_id or f"secuencia_{len(records) + 1}"),
                "description": current_desc,
                "type": kind,
                "length": len(seq),
                "gc_content": gc_content(seq, kind),
                "sequence": seq,
            })

    for raw in content.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            header = line[1:].strip()
            parts = header.split(None, 1)
            current_id = parts[0] if parts else ""
            current_desc = parts[1] if len(parts) > 1 else ""
            current_seq = []
        else:
            current_seq.append(line)

    flush()

    if not records:
        raise SequenceError("No se encontraron registros FASTA válidos.")
    return records