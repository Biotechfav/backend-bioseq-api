"""Enzimas de restricción y sus sitios de corte.

Las endonucleasas de restricción cortan el ADN en secuencias cortas y
específicas ("sitios de restricción"). Son la base de la biología
molecular: permiten clonar genes, verificar digestión de plásmidos y
hacer huellas de ADN. Aquí se localizan los sitios de corte de las
enzimas más comunes de laboratorio.
"""

from .sequences import SequenceError, validate_nucleic

# Reconocimiento expresado en la hebra sense (5'→3'); `cut` es el índice
# (0-based) dentro del sitio donde se corta esa hebra.
ENZYMES: dict[str, dict] = {
    "EcoRI":   {"site": "GAATTC",   "cut": 1, "overhang": "AATT"},
    "BamHI":   {"site": "GGATCC",   "cut": 1, "overhang": "GATC"},
    "HindIII": {"site": "AAGCTT",   "cut": 1, "overhang": "AGCT"},
    "NotI":    {"site": "GCGGCCGC", "cut": 2, "overhang": "GGCC"},
    "SalI":    {"site": "GTCGAC",   "cut": 1, "overhang": "TCGA"},
    "XbaI":    {"site": "TCTAGA",   "cut": 1, "overhang": "CTAG"},
    "EcoRV":   {"site": "GATATC",   "cut": 3, "overhang": "AT"},
    "PstI":    {"site": "CTGCAG",   "cut": 5, "overhang": "TGCA"},
    "SacI":    {"site": "GAGCTC",   "cut": 5, "overhang": "GCTC"},
}

# Sitio de corte de la hebra complementaria (para reportarlo también 5'→3').
def _reverse_strand_cut(site: str, cut: int) -> int:
    """Índice (0-based) del corte en la hebra complementaria 5'→3'."""
    return len(site) - cut


def restriction_sites(sequence: str, enzymes: list[str] | None = None) -> dict:
    """Localiza los sitios de restricción en una secuencia de ADN.

    Devuelve, por enzima, cada posición de corte (0-based) sobre la
    hebra sense y la complementaria. Si `enzymes` es None, busca todas
    las enzimas del catálogo; si la lista trae un nombre desconocido,
    se indica en `errors`.
    """
    seq = validate_nucleic(sequence, "dna")
    requested = list(enzymes) if enzymes else list(ENZYMES)
    unknown = [e for e in requested if e not in ENZYMES]
    if unknown:
        raise SequenceError(f"Enzimas no reconocidas: {', '.join(unknown)}. "
                            f"Disponibles: {', '.join(ENZYMES)}")

    result: dict[str, dict] = {}
    for name in requested:
        info = ENZYMES[name]
        site, cut = info["site"], info["cut"]
        cuts = []
        start = 0
        while True:
            pos = seq.find(site, start)
            if pos == -1:
                break
            cuts.append({
                "position": pos,
                "site": site,
                "sense_cut": pos + cut,
                "reverse_cut": pos + _reverse_strand_cut(site, cut),
            })
            start = pos + 1
        result[name] = {"site": site, "cuts": cuts}

    return result