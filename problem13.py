import sys
sys.path.insert(0, 'C:/Users/ADMIN/Desktop/Coder/PythonFEM/Functions')
from outputDisplacementReactions import *
import numpy as np
# Khai báo số liệu cần nhập
# E: Module đàn hồi
# I: Moment thứ 2 của khu vực
# L: Chiều dài của thanh
E = 210e6
A = 0.02
Iy = 10e-5
Iz = 20e-5
J = 5e-5
G = 84e6

# generation of coordinates and connectivities
nodeCoordinates = np.matrix([[0,0,0],[0,0,4],[4,0,4],[4,0,0],[0,5,0],[0,5,4],[4,5,4],[4,5,0]])
xx = np.zeros((1, np.size(nodeCoordinates,0)))
yy = np.zeros((1, np.size(nodeCoordinates,0)))
zz = np.zeros((1, np.size(nodeCoordinates,0)))
for i in range (0,len(nodeCoordinates)):
    xx[0, i]=nodeCoordinates[i,0]
    yy[0, i]=nodeCoordinates[i,1]
    zz[0, i]=nodeCoordinates[i, 2]
elementNodes = np.matrix ([[1,5],[2,6],[3,7],[4,8],[5,6],[6,7],[7,8],[8,5]])
numberNodes = np.size(nodeCoordinates, 0)
numberElements = np.size(elementNodes, 0)

# Cho cấu trúc
# displacements: vector chuyển vị
#  force: vector lực
# stiffness: ma trận độ cứng
# GDof: tổng số bậc tự do

GDof = 6*numberNodes
U=np.zeros((GDof,1))
force = np.zeros((GDof, 1))
stiffness = np.zeros((GDof, GDof))

# Vector lực
force[36] = -15

# Ma trận độ cứng
from formStiffness3Dframe import *
[stiffness,R,k,elementDof] = formStiffness3Dframe(GDof, numberElements,elementNodes, numberNodes, nodeCoordinates, E, A, Iz, Iy, G, J)
prescribedDof=[]
for i in range (1,25):
    prescribedDof.append([i])
prescribedDof=np.matrix(prescribedDof)
from solution import *
displacements = solution(GDof, prescribedDof, stiffness, force)

outputDisplacementReactions(displacements,stiffness,GDof,prescribedDof)