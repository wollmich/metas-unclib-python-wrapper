# Example Gum H1 - End-gauge calibration
# Michael Wollensack METAS - 12.02.2019 - 15.03.2019

print('Example Gum H1 - End-gauge calibration')
print('Begin')

from metas_unclib import *
use_linprop()

# Calibration of standard end gauge
l_s = ufloat(50.000623e-3, 25e-9, 1./18.)

# Measured difference between end gauges
# repeated observations
d1 = ufloat(215e-9, 5.8e-9, 1./24.)
# random effects of comporator
d2 = ufloat(0, 3.9e-9, 1./5.)
# systematic effects of comporator
d3 = ufloat(0, 6.7e-9, 1./8.)
d = d1 + d2 + d3
print('d       = %0.6e m' % get_value(d))
print('u(d)    = %0.6e m' % get_stdunc(d))
print('dof(d)  = %0.2f' % (1./get_idof(d)))

# Thermal expansion coefficient of standard end gauge (uniform)
alpha_s = ufloat(11.5e-6, 2e-6/np.sqrt(3))

# Temperature of test bed
# mean temperature of bed
theta_1 = ufloat(-0.1, 0.2)
# cyclic variation of temperature of room (arcsine)
theta_2 = ufloat(0, 0.5/np.sqrt(2))
theta = theta_1 + theta_2

# Difference in expansion coefficients of end gauges (uniform)
delta_alpha = ufloat(0, 1e-6/np.sqrt(3), 1./50.)

# Difference in temperatures of end gauges (uniform)
delta_theta = ufloat(0, 0.05/np.sqrt(3), 1./2.)

# Mathematical model 1
alpha = delta_alpha + alpha_s
theta_s = theta - delta_theta
l1 = (l_s * (1 + alpha_s * theta_s) + d) / (1 + alpha * theta)
print('Final result:')
print('l1      = %0.6e m' % get_value(l1))
print('u(l1)   = %0.6e m' % get_stdunc(l1))
print('dof(l1) = %0.2f' % (1./get_idof(l1)))

# Mathematical model 2
tmp1 = -l_s * delta_alpha * theta
tmp2 = -l_s * alpha_s * delta_theta
l2 = l_s + d + tmp1 + tmp2
print('Final result:')
print('l2      = %0.6e m' % get_value(l2))
print('u(l2)   = %0.6e m' % get_stdunc(l2))
print('dof(l2) = %0.2f' % (1./get_idof(l2)))

# Other
#print(get_correlation([l1, l2]))
#print(get_unc_component(l1, [l_s, d, alpha_s, theta, delta_alpha, delta_theta]))

print('End')
