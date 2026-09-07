"""Pruebas de integración de la API REST (FastAPI + SQLite)."""

import os
import tempfile

import pytest
from fastapi.testclient import TestClient

# Base de datos temporal para no ensuciar la de desarrollo.
os.environ["BIOSEQ_DB_PATH"] = os.path.join(tempfile.mkdtemp(), "test.db")

from app.main import app  # noqa: E402

client = TestClient(app)

GFP = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCA"


def test_health():
    r = client.get("/health")
    assert r.status_code == 200
    assert r.json()["status"] == "ok"


def test_root_info():
    r = client.get("/")
    assert r.status_code == 200
    assert "bioseq" not in r.json()  # solo comprobamos 200 y endpoints
    assert "endpoints" in r.json()


def test_composition_endpoint():
    r = client.post("/analyze/composition", json={"sequence": "ATGCGCTA"})
    assert r.status_code == 200
    data = r.json()
    assert data["length"] == 8
    assert data["gc_content"] == 50.0
    assert data["reverse_complement"] == "TAGCGCAT"


def test_translate_endpoint():
    r = client.post(
        "/analyze/translate",
        json={"sequence": GFP, "frame": 0},
    )
    assert r.status_code == 200
    assert r.json()["protein"].startswith("MVSKGEELFTG")


def test_translate_reverse_strand():
    r = client.post(
        "/analyze/translate",
        json={"sequence": "ATGGTGAGCAAG", "strand": "reverse", "frame": 0},
    )
    assert r.status_code == 200
    assert r.json()["protein"] == "LAHH"


def test_orfs_endpoint():
    r = client.post(
        "/analyze/orfs",
        json={"sequence": "ATGAAATAAATGCCCTAA", "min_aa": 1},
    )
    assert r.status_code == 200
    data = r.json()
    assert data["count"] == 2
    assert data["orfs"][0]["protein"] == "MK*"


def test_restriction_endpoint():
    r = client.post(
        "/analyze/restriction",
        json={"sequence": "GAATTCAATTGAATTC", "enzymes": ["EcoRI"]},
    )
    assert r.status_code == 200
    assert len(r.json()["EcoRI"]["cuts"]) == 2


def test_kmer_endpoint():
    r = client.post(
        "/analyze/kmer",
        json={"sequence": "AAAACCCC", "k": 2},
    )
    assert r.status_code == 200
    assert r.json()["counts"]["AA"] == 3


def test_align_endpoint():
    r = client.post(
        "/analyze/align",
        json={"seq_a": "GATTACA", "seq_b": "GATACA"},
    )
    assert r.status_code == 200
    data = r.json()
    assert data["score"] == 10
    assert data["identity_percent"] == 85.71


def test_fasta_endpoint():
    r = client.post("/fasta/parse", json={"content": ">a\nATGC\n>b\nGGGG"})
    assert r.status_code == 200
    assert r.json()["count"] == 2


def test_invalid_sequence_returns_422():
    r = client.post("/analyze/composition", json={"sequence": "ATZX"})
    assert r.status_code == 422


def test_unknown_endpoint_404():
    r = client.get("/nonexistent")
    assert r.status_code == 404


def test_history_recorded_and_queried():
    client.post("/analyze/composition", json={"sequence": "ATGC"})
    client.delete("/history")  # limpia antes de medir
    client.post("/analyze/composition", json={"sequence": "ATGC"})
    r = client.get("/history")
    data = r.json()
    assert data["total"] == 1
    assert data["items"][0]["tool"] == "composition"


def test_clear_history():
    client.post("/analyze/composition", json={"sequence": "ATGC"})
    r = client.delete("/history")
    assert r.status_code == 200
    assert client.get("/history").json()["total"] == 0