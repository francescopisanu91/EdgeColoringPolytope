import edgeRelaxation
import verticesFinder
import numpy as np 
import testBoxTDI
import testParallelBoxTDI
import testIDP

# BEFORE PROCEDING CHOSE boxTDI TEST FUN at the bottom of this page
# G = {'V':[0,1,2],'E':[(0,1),(1,2),(0,2)]} #K3 --- Not boxTDI TEST FUN 2 has IDP
# G = {'V':[0,1,2,3],'E':[(0,1),(1,2),(0,2),(0,3)]} #K3+antenna --- Not boxTDI TEST FUN 2 has IDP
G = {'V':[0,1,2,3,4,5,6],'E':[(0,1),(0,2),(0,3),(1,2),(1,4),(2,5),(3,4),(3,5),(4,5)]} #bar C6 --- Not BoxTDI TEST FUN 1
# G = {'V':[0,1,2,3,4,5],'E':[(0,1),(1,2),(2,3),(3,4),(4,5),(0,5)]} #C6 --- BoxTDI TEST FUN 2
# G = {'V':[0,1,2,3,4],'E':[(0,1),(1,2),(2,3),(0,3),(0,4)]} #Q --- Not BoxTDI TEST FUN 2 has IDP
# G = {'V':[0,1,2,3],'E':[(0,1),(1,2),(2,3),(0,3)]} #C4 --- BoxTDI TEST FUN 1&2 has IDP

G_edges = G['E']
n_edges = len(G['E'])
n_vertices = len(G['V'])
k = 3

edgeRelxationMatrix = edgeRelaxation.edgeRelaxationEdgeColoringMatrix(G_edges, n_vertices, k)
# print(edgeRelxationMatrix)
edgeRelaxationRHS = edgeRelaxation.edgeRelaxationEdgeColoringRHS(n_vertices, n_edges, k)
# print(edgeRelaxationRHS)

vertices = verticesFinder.polyVertices(edgeRelxationMatrix, edgeRelaxationRHS)
# for v in vertices:
#     print(v)

integrals, fractionals = verticesFinder.split_integral_fractional(vertices)

# for v in integrals:
#     print(v)

for v in fractionals:
    print(v)    
    # print(testIDP.find_minimal_integer_combination(2*v, integrals, max_coeff= k))
print(testIDP.hasIDP(2*np.asarray(fractionals), integrals, max_coeff= k))

# testBoxTDI.isPolyhedronBoxTDI1(edgeRelxationMatrix,edgeRelaxationRHS,vertices)
# testBoxTDI.isPolyhedronBoxTDI2(edgeRelxationMatrix,edgeRelaxationRHS,vertices)
# testParallelBoxTDI.isBoxTDI_parallel(edgeRelxationMatrix,edgeRelaxationRHS,vertices)