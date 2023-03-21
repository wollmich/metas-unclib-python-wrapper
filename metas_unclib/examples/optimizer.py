# Example Optimizer
# Michael Wollensack METAS - 13.02.2019 - 15.03.2019

print('Example Optimizer')
print('Begin')

from metas_unclib import *
use_linprop()

# Example Objective Function
def example_objective_function(x, p):
	print(get_value(x))
	f1 = p[0] * x[0] * x[0] + p[1] * x[0] + p[2]
	f2 = umath.exp(x[1]) - p[3]
	f = np.array([f1, f2])
	return f


# Start values for the variable optimization parameters
xstart = np.array([2.5, 2.0])
# Constant optimization parameters
p = np.array([ufloat(6.0, 2.0), ufloat(-4.0, 1.0), ufloat(-16.0, 1.0), ufloat(2.0, 1.0)])

# Analytic solution for example objective function (just for comparison)
x0 = np.array([[(-p[1] + umath.sqrt(p[1] * p[1] - 4.0 * p[0] * p[2])) / (2.0 * p[0])], [umath.log(p[3])]])

# Optimize using Metas.UncLib.Optimization (Levenberg-Marquardt)
x1 = unumlib.optimizer(lambda _x, _p: example_objective_function(_x, _p), xstart, p)

# Optimize using Metas.UncLib.Optimization (Trust-Region)
x2 = unumlib.optimizer(lambda _x, _p: example_objective_function(_x, _p), xstart, p, epsx=1e-6, algorithm=OptimizerAlgorithm.TrustRegion)

print('End')
