# Example Numerical Step
# Michael Wollensack METAS - 15.03.2019

print('Example Numerical Step')
print('Begin')

from metas_unclib import *
use_linprop()

# Example Numerical Function
def example_numerical_function(x):
	print(x)
	a = x[0]
	b = x[1]
	c = np.sqrt(a * a + b * b)
	d = a + b + c
	y = np.array([c, d])
	return y

# Uncertainty inputs
a = ufloat(3, 0.3)
b = ufloat(4, 0.4)
# Compute result using linear uncertainty propagation
c = umath.sqrt(a * a + b * b)
d = a + b + c
# Compute result using a numerical function Y = f(X)
temp = unumlib.numerical_step(example_numerical_function, [a, b], [0.001, 0.001])
c2 = temp[0]
d2 = temp[1]
# Compare the results
jacobi = get_jacobi2([c, c2, d, d2], [a, b])
print(jacobi)

print('End')
