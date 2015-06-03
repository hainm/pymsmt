from pymsmtmol.gauio import get_matrix_from_fchk
import numpy

"""
fcmatrix = get_matrix_from_fchk('H2O_opt.fchk', 'Cartesian Force Constants', 'Dipole Moment', 9)

fcmatrix = numpy.array([[round(j, 2) for j in i] for i in fcmatrix])

print fcmatrix


bfcmatrix = numpy.array([[float(0) for x in range(3)] for x in range(3)])

for i in range(0, 3):
  for j in range(0, 3):
    bfcmatrix[i][j] = -fcmatrix[3*(2-1)+i][3*(3-1)+j]

print bfcmatrix
"""

from numpy.linalg import eigvals, eig, norm
from mcpb.gene_final_frcmod_file import eigensort

a = numpy.array([[201, 0, 4], [0, 164, 0], [-1, 0, 231]])

at = numpy.transpose(a)

bfcmatrix = at

#Get the enigen vector and eigen value
eigval, eigvector = eig(bfcmatrix)
eigval, eigvector = eigensort(eigval, eigvector)

vec12 = [1.0/numpy.sqrt(3), 1.0/numpy.sqrt(3), 1.0/numpy.sqrt(3)]

fc1 = 0.0
fc2 = 0.0

for i in range(0, 3):
  ev = numpy.array([eigvector[0][i], eigvector[1][i], eigvector[2][i]])
  #if abs(numpy.dot(ev, vec12)) > 0.95:
  #If the direction of the force is almost the same with the
  #unit vector
  #  fc1 = eigval[i]
  fc2 = fc2 + eigval[i] * abs(numpy.dot(ev, vec12))

#if (fc1 > fc2):
#  fcfinal = fc1
#else:
fcfinal = fc2

print fcfinal

