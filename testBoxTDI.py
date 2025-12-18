import numpy as np
import itertools
import time
from math import comb

def reduceSize(A, tol=1e-9):
    """
    Rimuove colonne duplicate a meno del segno.
    Restituisce una matrice con colonne uniche (up to sign).
    """
    A = np.asarray(A, dtype=float)
    r, m = A.shape

    kept_cols = []
    seen = set()

    for j in range(m):
        col = A[:, j]

        # normalizzazione canonica (segno)
        # rendiamo la prima entry non-zero positiva
        if np.linalg.norm(col) < tol:
            key = tuple(np.zeros(r))
        else:
            for i in range(r):
                if abs(col[i]) > tol:
                    sign = 1 if col[i] > 0 else -1
                    break
            key = tuple(np.round(sign * col, 12))

        if key not in seen:
            seen.add(key)
            kept_cols.append(j)

    return A[:, kept_cols]

def isEquimodular(A, tol=1e-9):
    A = reduceSize(A, tol)
    A = np.asarray(A, dtype=float)
    r, m = A.shape

    if r > m:
        raise ValueError("La matrice deve avere r <= m")

    abs_dets = []

    # Tutte le scelte di r colonne
    for cols in itertools.combinations(range(m), r):
        sub = A[:, cols]
        det = np.linalg.det(sub)

        if abs(det) > tol:
            abs_dets.append(abs(det))

    # Se tutti i minori sono zero → equimodulare
    if len(abs_dets) == 0:
        return True

    # Tutti uguali (entro tolleranza)
    first = abs_dets[0]
    for d in abs_dets[1:]:
        if abs(d - first) > tol:
            return False

    return True

def rank(M, tol=1e-9):
    return np.linalg.matrix_rank(M, tol)

def print_progress(start_time, total, checked, bar_width=40):
        frac = checked / total
        filled = int(bar_width * frac)
        bar = "█" * filled + "-" * (bar_width - filled)
        elapsed = time.time() - start_time
        eta = elapsed * (total - checked) / checked if checked > 0 else 0
        print(
            f"\rProgress [{bar}] {100*frac:6.2f}% "
            f"({checked}/{total}) ETA {eta:7.1f}s",
            end="",
            flush=True
        )

def isPolyhedronBoxTDI(A, b, vertices, tol=1e-9, bar_width=40):
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float)
    V = [np.asarray(v, dtype=float) for v in vertices]

    n, m = A.shape
    at_most = min(n,m)

    # numero totale di sottoinsiemi di righe
    total = sum(comb(n, r) for r in range(1, at_most + 1))
    checked = 0

    tested_faces = []
    start_time = time.time()

    
    for r in range(1, at_most + 1):
        for rows in itertools.combinations(range(n), r):
            checked += 1
            if checked % 50 == 0 or checked == total:
                print_progress(start_time, total, checked, bar_width)

            A_sub = A[list(rows), :]
            b_sub = b[list(rows)]

            # (1) full row rank
            rank_A = np.linalg.matrix_rank(A_sub, tol)
            if  rank_A < r:
                continue
            
            # (2) vertici che saturano A_sub x = b_sub
            face_vertices = []
            face_indices = []

            for i, v in enumerate(V):
                if np.all(np.abs(A_sub @ v - b_sub) <= tol):
                    face_vertices.append(v)
                    face_indices.append(i)

            if not face_vertices:
                continue

            face_set = frozenset(face_indices)

            # skip facce già viste
            if face_set in tested_faces:
                continue

            tested_faces.append(face_set)

            # (3) almeno m - r + 1 aff. indipendenti
            needed = m - r + 1
            found = False

            for pts in itertools.combinations(face_vertices, needed):
                P = np.array(pts)
                D = P[1:] - P[0]
                if np.linalg.matrix_rank(D, tol) == needed - 1:
                    found = True
                    break

            if not found:
                continue

            # (4) equimodularità
            if not isEquimodular(A_sub):
                print("\nThe face given by")
                print("A_sub =")
                print(A_sub)
                print("b_sub =")
                print(b_sub)
                print("is not boxTDI")
                return False

    print("\nIl poliedro è boxTDI")
    return True

def isPolyhedronBoxTDI2(A, b, vertices, tol=1e-9, bar_width = 40):
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float)
    V = [np.asarray(v, dtype=float) for v in vertices]

    n, m = A.shape
    at_most = min(n,m)
    total = sum(comb(n, r) for r in range(1, at_most + 1))
    checked = 0
    start_time = time.time()

    for r in range(1, at_most + 1):
        for rows in itertools.combinations(range(n), r):
            checked += 1
            if checked % 50 == 0 or checked == total:
                print_progress(start_time, total, checked, bar_width)

            A_sub = A[list(rows), :]
            b_sub = b[list(rows)]

            # (1) full row rank
            if np.linalg.matrix_rank(A_sub, tol) < r:
                continue

            # (2) prima: test equimodularità
            if isEquimodular(A_sub):
                continue

            # (3) NON equimodulare → cerco punti che saturano
            face_vertices = []

            for v in V:
                if np.all(np.abs(A_sub @ v - b_sub) <= tol):
                    face_vertices.append(v)

            needed = m - r + 1
            if len(face_vertices) < needed:
                continue

            # controllo indipendenza affine
            found = False
            witness = None

            for pts in itertools.combinations(face_vertices, needed):
                P = np.array(pts)
                D = P[1:] - P[0]
                if np.linalg.matrix_rank(D, tol) == needed - 1:
                    found = True
                    witness = pts
                    break

            # (4) certificato di NON boxTDI
            if found:
                print("Il poliedro NON è boxTDI")
                print("A_sub =")
                print(A_sub)
                print("b_sub =")
                print(b_sub)
                print(f"Punti che soddisfano A_sub x = b_sub (almeno {needed} aff. indipendenti):")
                for p in witness:
                    print(p)
                return False

    print("Il poliedro è boxTDI")
    return True

