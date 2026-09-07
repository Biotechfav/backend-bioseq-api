"""Demo local: ejercita todos los endpoints de la BioSeq API.

Uso:  .venv/bin/python scripts/demo.py
Requiere el servidor corriendo (uvicorn) o usa TestClient.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fastapi.testclient import TestClient  # noqa: E402

from app.main import app  # noqa: E402

client = TestClient(app)

GFP = (
    "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAG"
    "CTGGACGGCGACGTAAACGGCCACAAGTTCAGC"
)

OLIGO = "ATGCTGACGTACGTCAGCTAGCTA"  # primer de ~23 pb, 50% GC

CASOS = [
    ("COMPOSICIÓN (bases, %GC, Tm, complementaria)", "/analyze/composition",
     {"sequence": OLIGO}),
    ("TRADUCCIÓN a proteína (GFP)", "/analyze/translate",
     {"sequence": GFP, "frame": 0}),
    ("ORF DE OPEN READING FRAMES", "/analyze/orfs",
     {"sequence": GFP + "TAA" * 20, "min_aa": 5}),
    ("ENZIMAS DE RESTRICCIÓN (EcoRI en el MCS)", "/analyze/restriction",
     {"sequence": "GAATTCATGGATCCGAATTC", "enzymes": ["EcoRI", "BamHI"]}),
    ("K-MERS (k=3)", "/analyze/kmer",
     {"sequence": "ATGATGATGCCC", "k": 3}),
    ("ALINEAMIENTO global (Needleman-Wunsch)", "/analyze/align",
     {"seq_a": "MVSKGEELFTGVVPILVELDGDV", "seq_b": "MVSKGEELFTGVVPILVELDGDV"}),
    ("ALINEAMIENTO local (Smith-Waterman)", "/analyze/align",
     {"seq_a": "ACGTGGC", "seq_b": "GGCACGT", "local": True}),
    ("PARSEO FASTA (ADN + ARN)", "/fasta/parse",
     {"content": ">glucosa-6-fosfato mutasa\naugcgu\n>desconocida\nacgtacgt"}),
]


def main() -> None:
    for titulo, path, payload in CASOS:
        print(f"\n=== {titulo} ===")
        r = client.post(path, json=payload)
        print("status:", r.status_code)
        data = r.json()
        if path == "/analyze/translate":
            print("proteína:", data["protein"])
        elif path == "/analyze/orfs":
            print(f"ORFs encontrados: {data['count']}")
            for o in data["orfs"]:
                print("  ", o["strand"], "f", o["frame"], o["start"],
                      "-", o["end"], "pb |", o["protein"])
        elif path == "/analyze/restriction":
            for enzima, info in data.items():
                print(f"  {enzima}: {len(info['cuts'])} corte(s)")
        elif path == "/analyze/composition":
            print("GC%:", data["gc_content"], "| Tm:", data["melting_temperature"],
                  "°C | revcomp:", data["reverse_complement"])
        elif path == "/analyze/kmer":
            top = "; ".join(f"{w}:{n}" for w, n in list(data["counts"].items())[:5])
            print("k:", data["k"], "| top:", top)
        elif path == "/analyze/align":
            print("score:", data["score"], "| identidad:",
                  data["identity_percent"], "%")
            print("  a:", data["aligned_a"])
            print("  b:", data["aligned_b"])
        elif path == "/fasta/parse":
            print("registros:", data["count"], "|",
                  " + ".join(f"{rec['id']}({rec['type']},{rec['length']}pb)"
                             for rec in data["records"]))

    print("\n=== HISTORIAL (SQLite) ===")
    r = client.get("/history")
    print("total:", r.json()["total"], "análisis:", r.json()["items"][:3])


if __name__ == "__main__":
    main()