import itertools
import numpy as np
from scipy.linalg import block_diag

def incidence_matrix(G_edges, n_vertices):
    m = len(G_edges)
    A = np.zeros((n_vertices, m), dtype=int)
    for j, (u, v) in enumerate(G_edges):
        A[u, j] = 1
        A[v, j] = 1
    return A

def A_k(G_edges, n_vertices, k):
    A = incidence_matrix(G_edges, n_vertices)
    return block_diag(*([A] * k))

def I_k(G_edges, k):
    m = len(G_edges)
    if k < 1:
        raise ValueError("k deve essere >= 1")
    I = np.eye(m)
    return np.hstack([I] * k)

def M_k(G_edges, n_vertices, k):
    A = A_k(G_edges, n_vertices, k)
    I = I_k(G_edges, k)
    return np.vstack([A, I])

def edgeRelaxationMatrix(G_edges, n_vertices, k):
    M = M_k(G_edges, n_vertices, k)
    m = len(G_edges)
    return np.vstack([M, -np.eye(k*m)])

def edgeRelaxationRHS(n_vertices,n_edges,k):
    return np.hstack([np.ones(k*n_vertices+n_edges),np.zeros(k*n_edges)])

