# Example Roots
# Michael Wollensack METAS - 15.02.2019 - 15.03.2019

print('Example Roots')
print('Begin')

from metas_unclib import *
use_linprop()

print('Roots of Quadratic Polynom')
# 3*x**2 - 2*x - 4 = 0

p2 = ufloat(3.0, 0.3, desc='p2')
p1 = ufloat(-2.0, 0.2, desc='p1')
p0 = ufloat(-4.0, 0.4, desc='p0')
p = np.array([p2, p1, p0])

r = unumlib.roots(p)
r_value = np.roots(get_value(p))
print(r)

print('Roots of Quartic Polynom')
# x**4 - 1 = 0
q4 = ufloat(1, 0.1, desc='q4')
q3 = ufloat(0, 0.1, desc='q3')
q2 = ufloat(0, 0.1, desc='q2')
q1 = ufloat(0, 0.1, desc='q1')
q0 = ufloat(-1, 0.1, desc='q0')
q = np.array([q4, q3, q2, q1, q0])

s = unumlib.roots(q)
s_value = np.roots(get_value(q))
print(s)

print('End')
