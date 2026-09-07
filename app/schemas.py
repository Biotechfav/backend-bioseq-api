"""Modelos de entrada/salida de la API (validación con Pydantic v2)."""

from pydantic import BaseModel, Field, field_validator

from .core import tools
from .core.sequences import SequenceError, validate_nucleic


class NucleicRequest(BaseModel):
    """Secuencia de ADN o ARN a analizar."""

    sequence: str = Field(..., description="Secuencia de bases (IUPAC)")
    kind: str = Field("dna", description="Tipo: 'dna' o 'rna'")

    @field_validator("sequence")
    @classmethod
    def _valid_seq(cls, value: str) -> str:
        try:
            return validate_nucleic(value, "dna")
        except SequenceError as exc:
            raise ValueError(str(exc)) from exc

    @field_validator("kind")
    @classmethod
    def _valid_kind(cls, value: str) -> str:
        value = value.lower()
        if value not in ("dna", "rna"):
            raise ValueError("kind debe ser 'dna' o 'rna'.")
        return value


class TranslationRequest(NucleicRequest):
    frame: int = Field(0, ge=0, le=2, description="Marco de lectura 0, 1 o 2")
    strand: str = Field("forward", description="'forward' o 'reverse'")

    @field_validator("strand")
    @classmethod
    def _valid_strand(cls, value: str) -> str:
        value = value.lower()
        if value not in ("forward", "reverse"):
            raise ValueError("strand debe ser 'forward' o 'reverse'.")
        return value


class TranslationResponse(BaseModel):
    frame: int
    length_aa: int
    protein: str
    has_stop_codon: bool


class ORFRequest(NucleicRequest):
    min_aa: int = Field(20, ge=1, description="Longitud mínima de la proteína")


class ORFResponse(BaseModel):
    strand: str
    frame: int
    start: int
    end: int
    length_bp: int
    length_aa: int
    sequence: str
    protein: str


class RestrictionRequest(NucleicRequest):
    enzymes: list[str] | None = Field(
        None,
        description="Enzimas a buscar (vacío/None = todas)",
        examples=[["EcoRI", "BamHI"]],
    )


class KmerRequest(NucleicRequest):
    k: int = Field(3, ge=1, le=12, description="Longitud de la palabra (k-mer)")


class AlignmentRequest(BaseModel):
    seq_a: str = Field(..., description="Primera secuencia")
    seq_b: str = Field(..., description="Segunda secuencia")
    local: bool = Field(False, description="True = Smith-Waterman (local)")
    match: int = Field(2, description="Puntaje de coincidencia")
    mismatch: int = Field(-1, description="Puntaje de no coincidencia")
    gap: int = Field(-2, description="Penalización por hueco (gap)")

    @field_validator("seq_a", "seq_b")
    @classmethod
    def _valid_seq(cls, value: str) -> str:
        try:
            return tools.canonical_protein(value)
        except SequenceError as exc:
            raise ValueError(str(exc)) from exc

    @field_validator("match", "mismatch", "gap")
    @classmethod
    def _score_fields(cls, value: int) -> int:
        return value


class FastaRequest(BaseModel):
    content: str = Field(
        ...,
        description="Texto FASTA (múltiples registros permitidos)",
        examples=[">seq1 descripcion\nACGTACGT"],
    )


class CompositionResponse(BaseModel):
    sequence: str
    length: int
    composition: dict[str, int]
    gc_content: float
    reverse_complement: str
    melting_temperature: float


class TranslationResponseItem(BaseModel):
    frame: int
    protein: str
    length_aa: int
    has_stop_codon: bool


# Reutilizable para las respuestas complejas.
ResponseOk = dict