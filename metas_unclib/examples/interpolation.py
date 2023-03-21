# Example Interpolation
# Michael Wollensack METAS - 13.02.2019 - 15.03.2019

print('Example Interpolation')
print('Begin')

from metas_unclib import *
use_linprop()

import scipy.special as special
import matplotlib.pyplot as plt

# Subfunction Residual
def residual(a, b):
	d = get_value(a[:]) - get_value(b[:])
	r = np.dot(d, d)
	return r

# Subfunction Test Interpolation
def test_interpolation(x, y, f):
	# Systematic uncertainty
	es = ufloat(0.0, 0.1)
	# Random uncertainty
	er = ufloatarray(np.zeros(x.size), 0.1 ** 2 * np.eye(x.size))
	# Function
	yb = y + es + er
	# New X
	xx = np.linspace(x[0], x[-1], 20 * (x.size - 1) + 1)
	# New Y (nominal function)
	yy = f(xx)
	# Linear interpolation
	yy_1a = unumlib.interpolation(x, yb, 1, xx)
	yy_1b = unumlib.interpolation2(x, yb, 1, xx)
	# Quadratic interpolation
	yy_2a = unumlib.interpolation(x, yb, 2, xx)
	yy_2b = unumlib.interpolation2(x, yb, 2, xx)
	# Cubic interpolation
	yy_3a = unumlib.interpolation(x, yb, 3, xx)
	yy_3b = unumlib.interpolation2(x, yb, 3, xx)
	# Spline interpolation (not-a-knot)
	yy_4a = unumlib.splineinterpolation(x, yb, xx, SplineBoundary.Not_a_Knot, 0.0, SplineBoundary.Not_a_Knot, 0.0)
	yy_4b = unumlib.splineinterpolation2(x, yb, xx, SplineBoundary.Not_a_Knot, 0.0, SplineBoundary.Not_a_Knot, 0.0)
	# Spline interpolation (natural spline)
	yy_5a = unumlib.splineinterpolation(x, yb, xx, SplineBoundary.Natural_Spline, 0.0, SplineBoundary.Natural_Spline, 0.0)
	yy_5b = unumlib.splineinterpolation2(x, yb, xx, SplineBoundary.Natural_Spline, 0.0, SplineBoundary.Natural_Spline, 0.0)
	# Spline interpolation (1st derivative)
	yy_6a = unumlib.splineinterpolation(x, yb, xx, SplineBoundary.First_Derivative, 0.0, SplineBoundary.First_Derivative, 0.0)
	yy_6b = unumlib.splineinterpolation2(x, yb, xx, SplineBoundary.First_Derivative, 0.0, SplineBoundary.First_Derivative, 0.0)
	# Spline interpolation (2nd derivative)
	# Same result as natural spline when 2nd derivatives values are set to zero
	yy_7a = unumlib.splineinterpolation(x, yb, xx, SplineBoundary.Second_Derivative, 0.0, SplineBoundary.Second_Derivative, 0.0)
	yy_7b = unumlib.splineinterpolation2(x, yb, xx, SplineBoundary.Second_Derivative, 0.0, SplineBoundary.Second_Derivative, 0.0)
	# Residuals (Interpolation vs Function)
	res_1 = residual(yy_1a, yy)
	res_2 = residual(yy_2a, yy)
	res_3 = residual(yy_3a, yy)
	res_4 = residual(yy_4a, yy)
	res_5 = residual(yy_5a, yy)
	res_6 = residual(yy_6a, yy)
	res_7 = residual(yy_7a, yy)
	# Names
	name_1 = 'Linear'
	name_2 = 'Quadratic'
	name_3 = 'Cubic'
	name_4 = 'Spline (not-a-knot)'
	name_5 = 'Spline (natural spline)'
	name_6 = 'Spline (1st derivative)'
	name_7 = 'Spline (2nd derivative)'
	# Figure
	subset = [0, 2, 3, 4];
	names = [name_1, name_2, name_3, name_4, name_5, name_6, name_7]
	yy_a = np.block([[yy_1a], [yy_2a], [yy_3a], [yy_4a], [yy_5a], [yy_6a], [yy_7a]]).T
	yy_b = np.block([[yy_1b], [yy_2b], [yy_3b], [yy_4b], [yy_5b], [yy_6b], [yy_7b]]).T
	
	fig = plt.figure(figsize=(12, 8))
	ax1 = plt.subplot(221)
	ax1.plot(x, get_value(yb), 'bo', xx, get_value(yy_a[:, subset]))
	ax1.set(title='Values of interpolation and spline')

	ax2 = plt.subplot(222)
	ax2.plot(x, get_value(yb), 'bo', xx, get_value(yy_b[:, subset]))
	ax2.set(title='Values of interpolation2 and spline2')

	ax3 = plt.subplot(223)
	ax3.plot(x, get_stdunc(yb), 'bo', xx, get_stdunc(yy_a[:, subset]))
	ax3.set(title='Uncertainties of interpolation and spline')

	ax4 = plt.subplot(224)
	ax4.plot(x, get_stdunc(yb), 'bo', xx, get_stdunc(yy_b[:, subset]))
	ax4.set(title='Uncertainties of interpolation2 and spline2')
	namses_subset = [names[i] for i in subset]
	namses_subset.insert(0, 'Points')
	ax4.legend(namses_subset, loc=4)

	plt.tight_layout()
	plt.show()

print('Interpolation of Sine Data')
x1 = np.array([0.0, 1.0, 2.5, 3.6, 5.0, 7.0, 8.1, 10.0])
y1 = np.sin(x1)
# Test interpolation
test_interpolation(x1, y1, lambda x: np.sin(x))

print('Interpolation of Distribution')
x2 = np.arange(-4., 5.)
y2 = np.array([0.0, 0.15, 1.12, 2.36, 2.36, 1.46, 0.49, 0.06, 0.0])
# Test interpolation
test_interpolation(x2, y2, lambda x: np.full(x.size, np.nan))

print('Interpolation of Data')
x3 = np.arange(-3., 4.)
y3 = np.array([-1.0, -1.0, - 1.0, 0.0, 1.0, 1.0, 1.0])
# Test interpolation
test_interpolation(x3, y3, lambda x: np.full(x.size, np.nan))

print('Interpolation of Oscillatory Sample Function')
x4 = np.arange(0., 26.)
y4 = special.j1(x4)
# Test interpolation
test_interpolation(x4, y4, lambda x: special.j1(x))

print('End')
