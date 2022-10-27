# Example Integrate
# Michael Wollensack METAS - 13.02.2019 - 15.03.2019

print('Example Integrate')
print('Begin')

from metas_unclib import *
use_linprop()

# Integral of sin(x) from 0 to pi/2
n1 = 5
n2 = 101
x1 = np.linspace(0, np.pi / 2., n1)
x2 = np.linspace(0, np.pi / 2., n2)

y1 = np.sin(x1)
y2 = np.sin(x2)

y1b = ufloatarray(y1, 1e-3 * np.eye(n1))
y2b = ufloatarray(y2, 1e-3 * np.eye(n2))

y3b = unumlib.interpolation(x1, y1b, 1, x2)
y4b = unumlib.interpolation(x1, y1b, 2, x2)
y5b = unumlib.interpolation(x1, y1b, 3, x2)
y6b = unumlib.splineinterpolation(x1, y1b, x2, SplineBoundary.Not_a_Knot, 0.0, SplineBoundary.Not_a_Knot, 0.0)

print('Integral of function')
r = -np.cos(np.pi / 2.) + np.cos(0.)
print(r)

print('Integral of values using trapz of n1 data points')
r1_0_end = np.trapz(y1, x1)
print(r1_0_end)

print('Integral (linear) of n1 data points')
r1_1 = unumlib.integrate(x1, y1b, 1)
r1_1_end = r1_1[-1]
d1_1_1_0_end = get_value(r1_1_end) - r1_0_end
print(r1_1_end)

print('Integral (quadratic) of n1 data points')
r1_2 = unumlib.integrate(x1, y1b, 2)
r1_2_end = r1_2[-1]
print(r1_2_end)

print('Integral (cubic) of n1 data points')
r1_3 = unumlib.integrate(x1, y1b, 3)
r1_3_end = r1_3[-1]
print(r1_3_end)

print('Integral (spline) of n1 data points')
r1_s = unumlib.splineintegrate(x1, y1b, SplineBoundary.Not_a_Knot, 0.0, SplineBoundary.Not_a_Knot, 0.0)
r1_s_end = r1_s[-1]
print(r1_s_end)

print('Integral of values using trapz of n2 data points')
r2_0_end = np.trapz(y2, x2)
print(r2_0_end)

print('Integral (linear) of n2 data points')
r2_1 = unumlib.integrate(x2, y2b, 1)
r2_1_end = r2_1[-1]
d2_1_2_0_end = get_value(r2_1_end) - r2_0_end
print(r2_1_end)

print('Integral (quadratic) of n2 data points')
r2_2 = unumlib.integrate(x2, y2b, 2)
r2_2_end = r2_2[-1]
print(r2_2_end)

print('Integral (cubic) of n2 data points')
r2_3 = unumlib.integrate(x2, y2b, 3)
r2_3_end = r2_3[-1]
print(r2_3_end)

print('Integral (spline) of n2 data points')
r2_s = unumlib.splineintegrate(x2, y2b, SplineBoundary.Not_a_Knot, 0.0, SplineBoundary.Not_a_Knot, 0.0)
r2_s_end = r2_s[-1]
print(r2_s_end)

print('Integral (linear) of interpolated (linear) data points')
r3_1 = unumlib.integrate(x2, y3b, 1)
r3_1_end = r3_1[-1]
d3_1_1_1_end = r3_1[-1] - r1_1[-1]
print(r3_1_end)

print('Integral (quadratic) of interpolated (quadratic) data points')
r4_2 = unumlib.integrate(x2, y4b, 2)
r4_2_end = r4_2[-1]
d4_2_1_2_end = r4_2[-1] - r1_2[-1]
print(r4_2_end)

print('Integral (cubic) of interpolated (cubic) data points')
r5_3 = unumlib.integrate(x2, y5b, 3)
r5_3_end = r5_3[-1]
d5_3_1_3_end = r5_3[-1] - r1_3[-1]
print(r5_3_end)

print('Integral (spline) of interpolated (spline) data points')
r6_s = unumlib.splineintegrate(x2, y6b, SplineBoundary.Not_a_Knot, 0.0, SplineBoundary.Not_a_Knot, 0.0)
r6_s_end = r6_s[-1]
d6_s_1_s_end = r6_s[-1] - r1_s[-1]
print(r6_s_end)

print('End')
