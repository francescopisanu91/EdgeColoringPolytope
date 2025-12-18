import edgeRelaxation
import verticesFinder
import numpy as np 
import testBoxTDI

# BEFORE PROCEDING CHOSE boxTDI TEST FUN at the bottom of this page
# G = {'V':[0,1,2,3,4,5,6],'E':[(0,1),(0,2),(0,3),(1,2),(1,4),(2,5),(3,4),(3,5),(4,5)]} #Not BoxTDI TEST FUN 1
# G = {'V':[0,1,2,3,4,5,6],'E':[(0,1),(1,2),(2,3),(3,4),(4,5),(0,5)]} #BoxTDI TEST FUN 2
G = {'V':[0,1,2,3,4],'E':[(0,1),(1,2),(2,3),(0,3),(0,4)]} #Not BoxTDI TEST FUN 2
# G = {'V':[0,1,2,3],'E':[(0,1),(1,2),(2,3),(0,3)]} #BoxTDI TEST FUN 1&2

G_edges = G['E']
n_edges = len(G['E'])
n_vertices = len(G['V'])
k = 2

edgeRelxationMatrix = edgeRelaxation.edgeRelaxationMatrix(G_edges, n_vertices, k)
# print(edgeRelxationMatrix)
edgeRelaxationRHS = edgeRelaxation.edgeRelaxationRHS(n_vertices, n_edges, k)
# print(edgeRelaxationRHS)

vertices = verticesFinder.polyVertices(edgeRelxationMatrix, edgeRelaxationRHS)
for v in vertices:
    print(v)

# testBoxTDI.isPolyhedronBoxTDI(edgeRelxationMatrix,edgeRelaxationRHS,vertices)
testBoxTDI.isPolyhedronBoxTDI2(edgeRelxationMatrix,edgeRelaxationRHS,vertices)
