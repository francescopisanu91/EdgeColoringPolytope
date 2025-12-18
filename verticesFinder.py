import cdd
import numpy as np

def polyVertices(A, b):
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)

    m, n = A.shape

    # Costruiamo la matrice H: [b | -A]
    H = np.hstack([b.reshape(-1, 1), -A])

    mat = cdd.matrix_from_array(H.tolist(), rep_type=cdd.RepType.INEQUALITY)
    poly = cdd.polyhedron_from_matrix(mat)
    gens = cdd.copy_generators(poly) 

    vertices = []
    for v in gens.array:
        if v[0] == 1:                 # if v[0]!=0 it is a ray not a vertex
            vertices.append(v[1:])
    return vertices


