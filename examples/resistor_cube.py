# Example Resistor Cube
# Michael Wollensack METAS - 12.02.2019 - 15.03.2019

print('Example Resistor Cube')
print('Begin')

from metas_unclib import *
use_linprop()

# Definition of input uncertainty objects R01 ... R12 and U
U   = 1.
R01 = ufloat(50.0, 0.1, desc='R01')
R02 = ufloat(50.0, 0.1, desc='R02')
R03 = ufloat(50.0, 0.1, desc='R03')
R04 = ufloat(50.0, 0.1, desc='R04')
R05 = ufloat(50.0, 0.1, desc='R05')
R06 = ufloat(50.0, 0.1, desc='R06')
R07 = ufloat(50.0, 0.1, desc='R07')
R08 = ufloat(50.0, 0.1, desc='R08')
R09 = ufloat(50.0, 0.1, desc='R09')
R10 = ufloat(50.0, 0.1, desc='R10')
R11 = ufloat(50.0, 0.1, desc='R11')
R12 = ufloat(50.0, 0.1, desc='R12')

# Kirchhoff's circuit laws --> linear equation system
Ux = np.array([ 0, 0, 0, 0, 0, 0, U, U, U, U, U, U ])

Rx = np.array([[  -1,   0,   0,   1,   0,   0,   0,   0,   1,   0,   0,   0 ],
               [   0,  -1,   0,   0,   1,   1,   0,   0,   0,   0,   0,   0 ],
               [   0,   0,  -1,   0,   0,   0,   1,   1,   0,   0,   0,   0 ],
               [   0,   0,   0,   1,   1,   0,   0,   0,   0,  -1,   0,   0 ],
               [   0,   0,   0,   0,   0,   1,   1,   0,   0,   0,  -1,   0 ],
               [   0,   0,   0,   0,   0,   0,   0,   1,   1,   0,   0,  -1 ],
               [ R01,   0,   0, R04,   0,   0,   0,   0,   0, R10,   0,   0 ],
               [ R01,   0,   0,   0,   0,   0,   0,   0, R09,   0,   0, R12 ],
               [   0, R02,   0,   0,   0, R06,   0,   0,   0,   0, R11,   0 ],
               [   0, R02,   0,   0, R05,   0,   0,   0,   0, R10,   0,   0 ],
               [   0,   0, R03,   0,   0,   0,   0, R08,   0,   0,   0, R12 ],
               [   0,   0, R03,   0,   0,   0, R07,   0,   0,   0, R11,   0 ]])

# Solve linear equation system
Ix = ulinalg.solve(Rx, Ux)

# Compute equivalent resistor of the cube
I = Ix[0] + Ix[1] + Ix[2]
R = U/I

print(R)

unc_budget(R, format='f6', name='R')

cv = get_covariance(Ix)
print(cv)

#Ux2 = ulinalg.dot(Rx, Ix)
#Ix2 = ulinalg.dot(ulinalg.inv(Rx), Ux.reshape(Ux.size, 1))

print('End')
