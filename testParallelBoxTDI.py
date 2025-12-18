import numpy as np
import itertools
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import threading
import time
from math import comb

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

def face_witness(A_sub, b_sub, vertices, tol=1e-9):
    """
    Ritorna:
    - None → NON è una faccia vera
    - tuple(sorted indices)) → faccia vera identificata dai punti
    """
    sats = []

    for i, v in enumerate(vertices):
        if np.all(np.abs(A_sub @ v - b_sub) <= tol):
            sats.append((i, v))

    r, m = A_sub.shape
    needed = m - r + 1

    if len(sats) < needed:
        return None

    for comb in itertools.combinations(sats, needed):
        pts = np.array([v for _, v in comb])
        D = pts[1:] - pts[0]
        if np.linalg.matrix_rank(D, tol) == needed - 1:
            return tuple(sorted(i for i, _ in comb))

    return None


def isBoxTDI_parallel(A, b, vertices, tol=1e-9, bar_width=40):
    A = np.asarray(A, float)
    b = np.asarray(b, float)
    vertices = [np.asarray(v, float) for v in vertices]

    n, m = A.shape
    at_most = min(n,m)
    seen_faces = set()

    total = sum(comb(n, r) for r in range(1, at_most + 1))
    checked = 0
    start_time = time.time()

    for r in range(1, at_most + 1):
        for rows in itertools.combinations(range(n), r):
            checked += 1
            if checked % 50 == 0 or checked == total:
                print_progress(start_time, total, checked, bar_width)

            A_sub = A[list(rows)]
            b_sub = b[list(rows)]

            # full row rank
            if np.linalg.matrix_rank(A_sub, tol) < r:
                continue

            stop_event = threading.Event()
            result = {}

            def equimodular_task():
                if stop_event.is_set():
                    return
                eq = isEquimodular(A_sub)
                result['equimodular'] = eq
                if eq:
                    stop_event.set()

            def face_task():
                if stop_event.is_set():
                    return
                face = face_witness(A_sub, b_sub, vertices, tol)

                if face is None:
                    result['face'] = False
                    stop_event.set()
                    return

                if face in seen_faces:
                    result['face'] = 'seen'
                    stop_event.set()
                    return

                result['face'] = face
                # NON fermiamo: serve sapere equimodularità

            with ThreadPoolExecutor(max_workers=10) as ex:
                f1 = ex.submit(equimodular_task)
                f2 = ex.submit(face_task)

                wait([f1, f2], return_when=FIRST_COMPLETED)

                # CASI DI SKIP
                if 'equimodular' in result and result['equimodular']:
                    continue

                if 'face' in result and result['face'] in (False, 'seen'):
                    continue

                # Serve il secondo risultato
                wait([f1, f2])

                eq = result.get('equimodular', True)
                face = result.get('face', None)

                if eq is False and isinstance(face, tuple):
                    print("❌ Il poliedro NON è boxTDI")
                    print("A_sub =")
                    print(A_sub)
                    print("b_sub =")
                    print(b_sub)
                    print("Vertici della faccia:")
                    for i in face:
                        print(vertices[i])
                    return False

                if isinstance(face, tuple):
                    seen_faces.add(face)

    print("✅ Il poliedro è boxTDI")
    return True
