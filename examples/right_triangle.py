# Example Right Triangle
# Michael Wollensack METAS - 11.02.2019 - 15.03.2019

print('Example Right Triangle')
print('Begin')

from metas_unclib import *
use_linprop()

a = ufloat(3.0, 0.3, desc='a')
b = ufloat(4.0, 0.4, desc='b')
c = (a * a + b * b).sqrt()
p = a + b + c
A = a * b / 2

f = ucomplex(a, b)
g = a + b * 1j
c2 = f.abs()
fangldedeg = f.angle(True)

print(c)
print(p)
print(A)

print(get_jacobi(c))

print(get_jacobi2([c, p, A], [a, b]))
print(get_unc_component([c, p, A], [a, b]))
print(get_correlation([c, p, A]))
print(get_coverage_interval([c, p, A]))
print(get_value([a, b, c, f]))
print(get_moment(f, 2))

print(f)
print(g)
print(c2)

unc_budget(c, format='f4', name='c')

unc_budget([c, p], format='f4', infos=['c', 'p'])

print('End')
