# Example Gum H3 - Calibration of a thermometer
# Michael Wollensack METAS - 13.02.2019 - 15.03.2019

print('Example Gum H3 - Calibration of a thermometer')
print('Begin')

from metas_unclib import *
use_linprop()

import matplotlib.pyplot as plt

# Data used to obtain a linear calibration curve for a thermometer
t0 = 20.
# Thermometer reading
tk = np.array([21.521, 22.012, 22.512, 23.003, 23.507, 23.999, 24.513, 25.002, 25.503, 26.010, 26.511])
# Observed correction
bk = np.array([-0.171, -0.169, -0.166, -0.159, -0.164, -0.165, -0.156, -0.157, -0.159, -0.161, -0.160])


# Least-square fitting (see numpy help polyfit)
p, cv = np.polyfit((tk - t0), bk, 1, cov=True)

tmp = ufloatarray(p, cv)
y1 = tmp[1]
y2 = tmp[0]
corr = get_correlation([y1, y2])

print('y1         = %10.6f °C' % get_value(y1))
print('u(y1)      = %10.6f °C' % get_stdunc(y1))
print('y2         = %10.6f' % get_value(y2))
print('u(y2)      = %10.6f' % get_stdunc(y2))
print('r(y1,y2)   = %10.6f' % corr[0, 1])

# Predicted correction and difference between observed and
b_tk = y1 + y2 * (tk - t0)
delta_b = bk - b_tk

# Uncertainty of predicted value (Example t = 30°C)
b_30 = y1 + y2 * (30. - t0)
print('b(30°C)    = %10.6f °C' % get_value(b_30))
print('u(b(30°C)) = %10.6f °C' % get_stdunc(b_30))

# Uncertainty of predicted value (Plot t = 20...30°C)
t = np.linspace(20., 30., 1001)
b = y1 + y2 * (t - t0)
value = get_value(b)
unc = get_stdunc(b)
v1 = value + unc
v2 = value - unc

fig, ax = plt.subplots(figsize=(12, 8))
ax.fill_between(t, v1, v2, facecolor='magenta', alpha=0.5)
ax.plot(t, value, 'b-', tk, bk, 'bo')

ax.set(xlabel='t (°C)', ylabel='b (°C)',
       title='Example Gum H3 - Calibration of a thermometer')
ax.grid()
plt.tight_layout()
plt.show()

print('End')
