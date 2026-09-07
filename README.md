# BioSeq API · Backend de Bioinformática

Proyecto **5 del portafolio backend** para **Practicante Backend — TEINOR S.A.C.**

API REST en **FastAPI** que procesa secuencias biológicas (ADN, ARN y
proteínas) resolviendo problemas reales de un laboratorio de
biotecnología con **algoritmos implementados en Python puro** — sin
librerías externas de biología: cada fórmula se programó a mano.

> Un backend que "entiende de biología": análisis de secuencias,
> traducción génica, búsqueda de genes (ORFs), enzimas de restricción,
> alineamiento por programación dinámica y manejo del formato FASTA.

## ¿Por qué es innovador?
No es otro CRUD: es un **servicio computacional de bioinformática**
(aplicado directo a la carrera de biotecnología del autor) con
algoritmos clásicos de la disciplina:

| Herramienta | Algoritmo / biología | Uso en laboratorio |
|-------------|----------------------|--------------------|
| %GC y Tm | Marmur-Doty | Diseño de primers de PCR |
| Traducción | Código genético estándar | Obtener proteína de un ORF |
| Detección de ORFs | Escaneo ATG→stop en 6 fases | Encontrar genes en ADN nuevo |
| Enzimas de restricción | Sitios de EcoRI, BamHI… | Planificar clonación |
| k-mers | Conteo de palabras | Huella genómica, identificación |
| Alineamiento | Needleman‑Wunsch / Smith‑Waterman | Comparar secuencias, homologías |
| FASTA | Parseo de registros | Formato NCBI/GenBank |

## Stack
- **FastAPI** + **Pydantic v2** (validación y schemas)
- **Uvicorn** (servidor ASGI)
- **SQLite** para el historial de análisis (sin servidor extra)
- **pytest** (37 pruebas: 22 de algoritmos + 15 de API)

## Estructura
```
bioseq-api/
├── app/
│   ├── main.py                  # App FastAPI + Swagger (/docs)
│   ├── schemas.py               # Modelos de entrada/salida (Pydantic)
│   ├── db.py                    # Historial SQLite
│   ├── core/                    # ← Ciencia (algoritmos puros)
│   │   ├── sequences.py         #   validación IUPAC, composición, %GC, Tm
│   │   ├── tools.py             #   traducción, ORFs, k-mers
│   │   ├── restriction.py       #   enzimas de restricción
│   │   ├── alignment.py         #   Needleman-Wunsch / Smith-Waterman
│   │   └── fasta.py             #   parseo FASTA
│   └── routers/
│       ├── analysis.py          # POST /analyze/*
│       └── util.py              # /fasta/parse, /history
├── tests/                       # pytest
├── scripts/demo.py              # Demo de todos los endpoints
└── requirements.txt
```

## Cómo ejecutar
```bash
python -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/uvicorn app.main:app --port 8000
```
- Documentación interactiva (Swagger): http://localhost:8000/docs
- JSON de la API: http://localhost:8000/openapi.json

## Endpoints (todos registran en /history)

| Método | Ruta | Descripción |
|--------|------|-------------|
| `GET`  | `/health` | Estado del servicio |
| `POST` | `/analyze/composition` | Bases, %GC, Tm, hebra complementaria |
| `POST` | `/analyze/translate` | ADN/ARN → proteína |
| `POST` | `/analyze/orfs` | Marcos de lectura abiertos (ATG→stop) |
| `POST` | `/analyze/restriction` | Cortes de enzimas de restricción |
| `POST` | `/analyze/kmer` | Frecuencia de k-mers |
| `POST` | `/analyze/align` | Alineamiento global o local |
| `POST` | `/fasta/parse` | Parseo de texto FASTA |
| `GET` / `DELETE` | `/history` | Ver / borrar historial |

### Ejemplo (cURL)
```bash
curl -s http://localhost:8000/analyze/composition \
  -H 'Content-Type: application/json' \
  -d '{"sequence":"ATGCTGACGTACGTCAGCTAGCTA"}'
```
```json
{
  "sequence": "ATGCTGACGTACGTCAGCTAGCTA",
  "length": 24,
  "composition": {"A": 6, "C": 6, "G": 6, "T": 6},
  "gc_content": 50.0,
  "reverse_complement": "TAGCTAGCTGACGTACGTCAGCAT",
  "melting_temperature": 57.38
}
```

```bash
# Alineamiento global de dos oligos
curl -s http://localhost:8000/analyze/align \
  -d '{"seq_a":"MVSKGEELFTGVVPILVELDGDV","seq_b":"MVSKGEELFTGVVPILVELDGDV"}'
```

## Validación y manejo de errores
- Validación estricta **IUPAC** (símbolos inválidos → `422`).
- Límites documentados (alineamiento ≤ 2000 pb, k ≤ 12).
- Tm fiable solo para secuencias ≥ 10 bases (devuelve `null` si no).
- Rutas desconocidas → `404`; errores de análisis → `422` con mensaje claro.

## Calidad
```bash
.venv/bin/python -m pytest tests/ -q   # 37 passed
.venv/bin/python scripts/demo.py        # demo de todos los endpoints
```

## Autor
Proyecto del portafolio de **Practicante Backend TEINOR S.A.C.** —
[GitHub: Biotechfav](https://github.com/Biotechfav)