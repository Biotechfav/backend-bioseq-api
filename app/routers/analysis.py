"""Routers de análisis bioinformático."""

from typing import Annotated

from fastapi import APIRouter, Depends, HTTPException

from .. import db
from ..core import alignment, restriction as restriction_core, sequences, tools
from ..schemas import (
    AlignmentRequest,
    KmerRequest,
    NucleicRequest,
    ORFRequest,
    RestrictionRequest,
    TranslationRequest,
)


router = APIRouter(prefix="/analyze", tags=["Análisis de secuencias"])


def _tool_error_handler(exc: Exception) -> HTTPException:
    if isinstance(exc, (sequences.SequenceError, ValueError)):
        return HTTPException(status_code=422, detail=str(exc))
    return HTTPException(status_code=500, detail=str(exc))


@router.post(
    "/composition",
    summary="Composición de bases, %GC, Tm y hebra complementaria",
    description=(
        "Analiza una secuencia de ADN/ARN: conteo de bases, contenido "
        "GC (estabilidad térmica), temperatura de melting estimada y "
        "cadena complementaria reversa."
    ),
)
def composition(body: NucleicRequest):
    try:
        seq = sequences.validate_nucleic(body.sequence, body.kind)
        result = {
            "sequence": seq,
            "length": len(seq),
            "composition": sequences.composition(seq, body.kind),
            "gc_content": sequences.gc_content(seq, body.kind),
            "reverse_complement": sequences.reverse_complement(seq, body.kind),
            "melting_temperature": sequences.melting_temperature(seq, body.kind),
        }
    except Exception as exc:  # noqa: BLE001
        raise _tool_error_handler(exc) from exc

    db.record("composition", len(seq))
    return result


@router.post(
    "/translate",
    summary="Traducir ADN/ARN a proteína",
    description=(
        "Traduce una secuencia a su proteína usando el código genético "
        "estándar. Indica el marco de lectura y si terminó en codón de "
        "parada (stop)."
    ),
)
def translate(body: TranslationRequest):
    try:
        seq = sequences.validate_nucleic(body.sequence, body.kind)
        if body.strand == "reverse":
            seq = sequences.reverse_complement(seq, body.kind)
        result = tools.translate(seq, body.kind, body.frame)
    except Exception as exc:  # noqa: BLE001
        raise _tool_error_handler(exc) from exc

    db.record("translate", len(seq), result["protein"][:40])
    return result


@router.post(
    "/orfs",
    summary="Detectar marcos de lectura abiertos (ORFs)",
    description=(
        "Busca ORFs en ambas hebras y en los 3 marcos de lectura: desde "
        "el codón de inicio ATG hasta el primer codón de parada. Es el "
        "primer paso para identificar una proteína codificada."
    ),
)
def orfs(body: ORFRequest):
    try:
        seq = sequences.validate_nucleic(body.sequence, body.kind)
        result = tools.find_orfs(seq, body.kind, body.min_aa)
    except Exception as exc:  # noqa: BLE001
        raise _tool_error_handler(exc) from exc

    db.record("orfs", len(seq), f"{len(result)} ORF(s)")
    return {"count": len(result), "orfs": result}


@router.post(
    "/restriction",
    summary="Sitios de corte de enzimas de restricción",
    description=(
        "Localiza los sitios de reconocimiento de enzimas de restricción "
        "comunes (EcoRI, BamHI, HindIII, NotI...) y la posición exacta "
        "de corte en ambas hebras."
    ),
)
def restriction(body: RestrictionRequest):
    try:
        seq = sequences.validate_nucleic(body.sequence, "dna")
        result = restriction_core.restriction_sites(seq, body.enzymes)
    except Exception as exc:  # noqa: BLE001
        raise _tool_error_handler(exc) from exc

    db.record("restriction", len(seq))
    return result


@router.post(
    "/kmer",
    summary="Frecuencia de k-mers",
    description=(
        "Calcula la frecuencia de cada palabra de k bases (k-mer) en la "
        "secuencia. Útil para huella genómica e identificación de especies."
    ),
)
def kmer(body: KmerRequest):
    try:
        seq = sequences.validate_nucleic(body.sequence, body.kind)
        result = tools.kmer_table(seq, body.k, body.kind)
    except Exception as exc:  # noqa: BLE001
        raise _tool_error_handler(exc) from exc

    total = sum(result.values())
    db.record("kmer", len(seq), f"k={body.k}, {total} palabras")
    return {"k": body.k, "total_kmers": total, "counts": result}


@router.post(
    "/align",
    summary="Alinear dos secuencias (Needleman-Wunsch / Smith-Waterman)",
    description=(
        "Alineamiento por programación dinámica: global "
        "(Needleman-Wunsch) o local (Smith-Waterman). Devuelve las "
        "cadenas alineadas y el porcentaje de identidad. Límite: 2000 "
        "caracteres por secuencia."
    ),
)
def align(body: AlignmentRequest):
    try:
        result = alignment.align(
            body.seq_a,
            body.seq_b,
            local=body.local,
            match=body.match,
            mismatch=body.mismatch,
            gap=body.gap,
        )
    except Exception as exc:  # noqa: BLE001
        raise _tool_error_handler(exc) from exc

    db.record("align" if not body.local else "align_local",
              max(len(body.seq_a), len(body.seq_b)),
              f"score={result['score']}, identidad={result['identity_percent']}%")
    return result