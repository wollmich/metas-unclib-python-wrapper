# Example ODE Solver
# Michael Wollensack METAS - 13.10.2025

print('Example ODE Solver')
print('Begin')

from metas_unclib import *
use_linprop()

# Example ODE Function
def example_ode_function(y, x, p ,dy):
    dy[0] = p[0] * y[0]

# Uncertainty inputs
y = [ufloat(1)]
x = [ufloat(0, 0.1), ufloat(1, 0.1), ufloat(2, 0.1), ufloat(3, 0.1)]
p = [ufloat(-1.2, 0.01)]
# Compute result
yres = [umath.exp(p[0] * xi) for xi in x]
# Compute result using ODE solver
ytbl = unumlib.ode_solver(y, x, p, 1e-12, example_ode_function)
# Compare the results
print(yres)

print(ytbl)

print('End')
