import numpy as np

def dominates(u, v):
    """ritorna True se u <= v componente per componente"""
    return np.all(u <= v)

def filter_pass_proper_coloring(integral_vertices, n_edges):
    """
    Tiene solo i vettori 0/1 con somma esattamente n_edges
    """
    filtered = []

    for v in integral_vertices:
        v = np.asarray(v, dtype=int)
        if v.sum() == n_edges:
            filtered.append(v)

    return filtered

def isPersistentFromBelow(integral_vertices, fractional_vertices, tol=1e-9):
    # assicuriamoci che i vertici interi siano array 0/1
    V_int = [np.asarray(v, dtype=int) for v in integral_vertices]

    for f in fractional_vertices:
        f = np.asarray(f, dtype=float)

        # floor (dato che i valori sono 0, 1/2, 1)
        u = np.floor(f + tol).astype(int)

        found = False
        for v in V_int:
            if dominates(u, v):
                found = True
                break

        if not found:
            print("❌ Test fallito")
            print("Vertice frazionario:", f)
            print("Floor:", u)
            return False

    print("✅ Test di persistenza superato")
    return True

def isPersistentFromAbove(integral_vertices, fractional_vertices, tol=1e-9):
    # assicuriamoci che i vertici interi siano array 0/1
    V_int = [np.asarray(v, dtype=int) for v in integral_vertices]

    for f in fractional_vertices:
        f = np.asarray(f, dtype=float)

        # ceil (dato che i valori sono 0, 1/2, 1)
        u = np.ceil(f - tol).astype(int)

        found = False
        for v in V_int:
            if dominates(v, u):
                found = True
                break

        if not found:
            print("❌ Test fallito")
            print("Vertice frazionario:", f)
            print("Ceil:", u)
            return False

    print("✅ Test di persistenza superato")
    return True
