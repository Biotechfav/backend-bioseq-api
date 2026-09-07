"""BioSeq API — herramientas de bioinformática.

Bootstrapping de FastAPI: documentación configurada para que Swagger sea
limpio y con ejemplos.
"""

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from .routers import analysis, util

DESCRIPTION = """
**BioSeq API** — plataforma backend de *bioinformática* que procesa
secuencias de ADN, ARN y proteínas con algoritmos implementados en
Python puro (sin librerías externas de biología).

### Herramientas disponibles
- 🧬 **Composición**: conteo de bases, %GC, temperatura de melting y
  hebra complementaria.
- 🔁 **Traducción**: ADN/ARN → proteína (código genético estándar).
- 🔎 **ORFs**: detección de marcos de lectura abiertos en ambas hebras.
- ✂️ **Enzimas de restricción**: sitios de corte (EcoRI, BamHI, …).
- 🧮 **k-mers**: frecuencia de palabras cortas.
- 🧬 **Alineamiento**: Needleman-Wunsch (global) y Smith-Waterman (local).
- 📄 **FASTA**: parseo del formato estándar de GenBank/NCBI.

Toda operación se registra en el **historial** (SQLite) consultable en
`GET /history`.
"""

app = FastAPI(
    title="BioSeq API · Backend de Bioinformática",
    version="1.0.0",
    description=DESCRIPTION,
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(analysis.router)
app.include_router(util.router)


@app.get("/", tags=["Info"], summary="Bienvenida y endpoints")
def root():
    return {
        "service": "BioSeq API",
        "version": "1.0.0",
        "swagger": "/docs",
        "health": "/health",
        "endpoints": {
            "POST /analyze/composition": "bases, %GC, Tm, complementaria",
            "POST /analyze/translate": "ADN/ARN → proteína",
            "POST /analyze/orfs": "marcos de lectura abiertos",
            "POST /analyze/restriction": "sitios de enzimas de restricción",
            "POST /analyze/kmer": "frecuencia de k-mers",
            "POST /analyze/align": "alineamiento NW / SW",
            "POST /fasta/parse": "parseo FASTA",
            "GET/DELETE /history": "historial de análisis",
        },
    }


@app.get("/health", tags=["Info"], summary="Estado del servicio")
def health():
    return {"status": "ok", "database": "sqlite"}