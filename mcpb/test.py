from numpy import average, array, dot, cross
from numpy.linalg import eigvals, eig, norm
from pymsmtmol.cal import calc_bond, calc_angle

crd1 = [-1.257875, 0.494274, 0.000000]
#crd2 = [-0.290487, 0.543584, 0.000000]
crd2 = [-1.534314, 1.422634, 0.000000]

disbohr = calc_bond(crd1, crd2) #unit is angstrom
vec12 = array(crd2) - array(crd1) #vec12 is vec2 - vec1
vec12 = [i/(disbohr) for i in vec12]
vec12 = array(vec12)

#bfcmatrix = [[0.382900e-04, 0.000000e+00, 0.000000e+00],
#           [0.000000e+00, -0.324222e+00, 0.186994e+00],
#           [0.000000e+00, 0.254875e+00, -0.226443e+00]]

#bfcmatrix = [[-0.511315, 0.000179, 0.000000],
#           [-0.060177, -0.044070, 0.000000],
#           [0.000000, 0.000000, 0.000708]]

bfcmatrix = [[-0.078101, 0.093541, 0.000000],
           [0.153187, -0.479088, 0.000000],
           [0.000000, 0.000000, -0.000116]]

bfcmatrix = array(bfcmatrix)
bfcmatrix = bfcmatrix * -1.0

eigval, eigvector = eig(bfcmatrix)
fc = 0.0
for i in range(0, 3):
  ev = eigvector[:,i]
  fc = fc + eigval[i] * abs(dot(ev, vec12))
fcfinal = fc * 2240.87 * 0.5
#2240.87 is Hatree/(Bohr^2) to kcal/(mol*angstrom^2)
#Times 0.5 factor since AMBER use k(r-r0)^2 but not 1/2*k*(r-r0)^2
#bfconst.append(fcfinal)
print fcfinal

