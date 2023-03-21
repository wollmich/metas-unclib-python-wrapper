# Example Eigenvalue Problem
# Michael Wollensack METAS - 15.02.2019 - 15.03.2019

print('Example Eigenvalue Problem')
print('Begin')

from metas_unclib import *
use_linprop()

# Symmetric Eigenvalue problem
s = ufloat(0.0, 0.1)
e11 = ufloat(1.0, 0.1)
e21 = ufloat(2.0, 0.1)
e22 = ufloat(5.0, 0.1)

m1 = np.array([[e11, e21], [e21, e22]]) + s

v1, d1 = ulinalg.eig(m1)
d1_value, v1_values = np.linalg.eig(get_value(m1))
# Check
check1 = ulinalg.dot(m1, v1) - ulinalg.dot(v1, np.diag(d1))


# Non-symmetric Eigenvalue Problem
m2 = np.array([[e11, 0.1 * e21], [e21, e22]]) + s

v2, d2 = ulinalg.eig(m2)
d2_value, v2_values = np.linalg.eig(get_value(m2))
# Check
check2 = ulinalg.dot(m2, v2) - ulinalg.dot(v2, np.diag(d2))


# Linear Eigenvalue Problem
a0 = np.array([[e11, e21], [e21, e22]]) + s
a1 = np.array([[1.0, 2.0], [0.0, 3.0]]) + s

v3, d3 = ulinalg.eig(a0, a1)
# Check
check3 = ulinalg.dot(a0, v3) + ulinalg.dot(a1, ulinalg.dot(v3, np.diag(d3)))


# Quadratic Eigenvalue Problem
a2 = np.array([[4.0, 0.0], [5.0, 6.0]]) + s

v4, d4 = ulinalg.eig(a0, a1, a2)
# Check
check4 = ulinalg.dot(a0, v4) + ulinalg.dot(a1, ulinalg.dot(v4, np.diag(d4))) + ulinalg.dot(a2, ulinalg.dot(v4, np.diag(d4**2)))


# Cubic Eigenvalue Problem
a3 = np.array([[7.0, 8.0], [9.0, 0.0]]) + s

v5, d5 = ulinalg.eig(a0, a1, a2, a3)
# Check
check5 = ulinalg.dot(a0, v5) \
	+ ulinalg.dot(a1, ulinalg.dot(v5, np.diag(d5))) \
	+ ulinalg.dot(a2, ulinalg.dot(v5, np.diag(d5**2))) \
    + ulinalg.dot(a3, ulinalg.dot(v5, np.diag(d5**3)))

# Over-determined Linear Eigenvalue Problem
b0 = np.array([[e11, e21], [e21, e22], [2.0, 3.0]]) + s
b1 = np.array([[4.0, 5.0], [6.0, 7.0], [8.0, 9.0]]) + s

v6, d6 = ulinalg.eig(b0, b1)
# Check
check6 = ulinalg.dot(ulinalg.dot(b0.T, b0), v6) \
    + ulinalg.dot(ulinalg.dot(b0.T, b1) + ulinalg.dot(b1.T, b0), ulinalg.dot(v6, np.diag(d6))) \
    + ulinalg.dot(ulinalg.dot(b1.T, b1), ulinalg.dot(v6, np.diag(d6**2)))

# Over-determined Quadratic Eigenvalue Problem
b2 = np.array([[10.0, 11.0], [12.0, 13.0], [14.0, 15.0]]) + s

v7, d7 = ulinalg.eig(b0, b1, b2)
# Check
check7 = ulinalg.dot(ulinalg.dot(b0.T, b0), v7) \
	+ ulinalg.dot(ulinalg.dot(b0.T, b1) + ulinalg.dot(b1.T, b0), ulinalg.dot(v7, np.diag(d7))) \
	+ ulinalg.dot(ulinalg.dot(b0.T, b2) + ulinalg.dot(b1.T, b1) + ulinalg.dot(b2.T, b0), ulinalg.dot(v7, np.diag(d7**2))) \
	+ ulinalg.dot(ulinalg.dot(b1.T, b2) + ulinalg.dot(b2.T, b1), ulinalg.dot(v7, np.diag(d7**3))) \
	+ ulinalg.dot(ulinalg.dot(b2.T, b2), ulinalg.dot(v7, np.diag(d7**4)))

print('End')
