"""Alineamiento de secuencias (programación dinámica).

Implementación didáctica de dos algoritmos clásicos de bioinformática:
  * Needleman-Wunsch (alineamiento global),
  * Smith-Waterman (alineamiento local).

Permiten comparar dos secuencias biológicas y visualizar su
similitud — base para buscar homologías, primers o mutaciones.
La complejidad es O(n·m) en tiempo y memoria, por eso el límite de
tamaño está fijado en 2000 caracteres por secuencia.
"""

from .sequences import normalize

MAX_ALIGNMENT_LENGTH = 2000


def _traceback(direction: list[list[str]], i: int, j: int,
               a: list[str], b: list[str]) -> tuple[str, str, int, int]:
    """Reconstruye el alineamiento desde la dirección de cada celda."""
    aligned_a: list[str] = []
    aligned_b: list[str] = []
    while i > 0 or j > 0:
        d = direction[i][j]
        if d == "diag":
            aligned_a.append(a[i - 1])
            aligned_b.append(b[j - 1])
            i -= 1
            j -= 1
        elif d == "up":
            aligned_a.append(a[i - 1])
            aligned_b.append("-")
            i -= 1
        elif d == "left":
            aligned_a.append("-")
            aligned_b.append(b[j - 1])
            j -= 1
        else:
            break  # Smith-Waterman: terminó el segmento local
    return (
        "".join(reversed(aligned_a)),
        "".join(reversed(aligned_b)),
        i,
        j,
    )


def align(seq_a: str, seq_b: str, local: bool = False,
          match: int = 2, mismatch: int = -1, gap: int = -2) -> dict:
    """Alinea dos secuencias (global o local) y reporta el puntaje.

    Parámetros: match/mismatch/gap son los puntajes del modelo de
    sustitución. Devuelve las cadenas alineadas, la identidad (%) y
    las posiciones de inicio en cada secuencia original.
    """
    a = normalize(seq_a)
    b = normalize(seq_b)

    if len(a) > MAX_ALIGNMENT_LENGTH or len(b) > MAX_ALIGNMENT_LENGTH:
        raise ValueError(
            f"Longitud máxima de alineamiento: {MAX_ALIGNMENT_LENGTH} "
            "caracteres por secuencia."
        )

    n, m = len(a), len(b)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    direction = [[""] * (m + 1) for _ in range(n + 1)]

    if not local:
        for i in range(1, n + 1):
            score[i][0] = score[i - 1][0] + gap
            direction[i][0] = "up"
        for j in range(1, m + 1):
            score[0][j] = score[0][j - 1] + gap
            direction[0][j] = "left"

    i_max, j_max, best = 0, 0, 0

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            sub = score[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mismatch)
            up = score[i - 1][j] + gap
            left = score[i][j - 1] + gap

            if local:
                cur = max(0, sub, up, left)
            else:
                cur = max(sub, up, left)

            if cur == sub:
                direction[i][j] = "diag"
            elif cur == up:
                direction[i][j] = "up"
            elif cur == left:
                direction[i][j] = "left"
            else:
                direction[i][j] = "stop"

            score[i][j] = cur
            if local and cur > best:
                best = cur
                i_max, j_max = i, j

    if not local:
        best = score[n][m]

    end_i = i_max if local else n
    end_j = j_max if local else m

    aligned_a, aligned_b, start_a, start_b = _traceback(
        direction, end_i, end_j,
        [c for c in a], [c for c in b],
    )

    pairs = sum(
        1 for x, y in zip(aligned_a, aligned_b) if x == y and x != "-"
    )
    identity = round(pairs / max(len(aligned_a), 1) * 100, 2)

    return {
        "mode": "local" if local else "global",
        "score": best,
        "identity_percent": identity,
        "aligned_a": aligned_a,
        "aligned_b": aligned_b,
        "start_a": start_a,
        "start_b": start_b,
    }