# Example Gum H2 - Simultaneous resistance and reactance measurement
# Michael Wollensack METAS - 12.02.2019 - 15.03.2019

print('Example Gum H2 - Simultaneous resistance and reactance measurement')
print('Begin')

from metas_unclib import *
use_linprop()

# Definition of the inputs

meas = np.array([[5.007, 19.663e-3, 1.0456],
                 [4.994, 19.639e-3, 1.0438],
                 [5.005, 19.640e-3, 1.0468],
                 [4.990, 19.685e-3, 1.0428],
		         [4.999, 19.678e-3, 1.0433]])

input_values = np.mean(meas, axis=0)
input_covar = np.cov(meas, rowvar=False) / meas.shape[0]

inputs = ufloatarray(input_values, input_covar)

v = inputs[0]
i = inputs[1]
phi = inputs[2]

# Compute the outputs

r = v / i * umath.cos(phi)
x = v / i * umath.sin(phi)
z = v / i

outputs = np.array([r, x, z])

print('Output values:')
output_values = get_value(outputs)
print(output_values)
print('Output standard uncertainties:')
output_stdunc = get_stdunc(outputs)
print(output_stdunc)
print('Output correlation:')
output_corr = get_correlation(outputs)
print(output_corr)

print('End')
