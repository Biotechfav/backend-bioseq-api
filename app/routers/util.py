"""FASTA e historial."""

from fastapi import APIRouter, HTTPException, Query

from .. import db
from ..core.fasta import parse_fasta
from ..core.sequences import SequenceError
from ..schemas import FastaRequest

router = APIRouter(tags=["Utilidades"])


@router.post(
    "/fasta/parse",
    summary="Analizar secuencias en formato FASTA",
    description=(
        "Recibe texto en formato FASTA (uno o más registros) y devuelve "
        "cada secuencia con su longitud y contenido GC. Es el formato "
        "estándar de bancos como NCBI/GenBank."
    ),
)
def parse(body: FastaRequest):
    try:
        records = parse_fasta(body.content)
    except SequenceError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc

    db.record("fasta", sum(r["length"] for r in records),
              f"{len(records)} registro(s)")
    return {"count": len(records), "records": records}


@router.get(
    "/history",
    summary="Historial de análisis ejecutados",
    description="Lista los últimos análisis realizados en la API.",
)
def history(limit: int = Query(20, ge=1, le=100)):
    return {"total": db.count(), "items": db.list_all(limit)}


@router.delete(
    "/history",
    summary="Borrar historial",
    description="Elimina todos los registros de análisis (no toca datos Biol.).",
)
def clear_history():
    return {"deleted": db.clear()}