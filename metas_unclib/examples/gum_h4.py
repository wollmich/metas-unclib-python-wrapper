# Example Gum H4 - Measurement of activity
# Michael Wollensack METAS - 12.02.2019 - 15.03.2019

print('Example Gum H4 - Measurement of activity')
print('Begin')

from metas_unclib import *
use_linprop()

# Counting data for determining the activity concentration
T0 = 60. # min
l = 1.25894e-4 # min^-1
ts = np.array([ 243.74,  984.53, 1723.87, 2463.17, 3217.56, 3956.83])
Cs = np.array([  15380,   14978,   14394,   13254,   12516,   11058])
tb = np.array([ 305.56, 1046.10, 1785.43, 2524.73, 3279.12, 4018.38])
Cb = np.array([   4054,    3922,    4200,    3830,    3956,    3980])
tx = np.array([ 367.37, 1107.66, 1846.99, 2586.28, 3340.68, 4079.94])
Cx = np.array([  41432,   38706,   35860,   32238,   29640,   26356])

Rx = (Cx - Cb) / T0 * np.exp(l * tx)
Rs = (Cs - Cb) / T0 * np.exp(l * ts)

data = np.block([[Rx], [Rs]]).T

# Analysis of data
input_values = np.mean(data, axis=0)
input_covar = np.cov(data, rowvar=False) / data.shape[0]

inputs = ufloatarray(input_values, input_covar)

Rx_m = inputs[0]
Rs_m = inputs[1]
corr = get_correlation([Rx_m, Rs_m])

print('Rx         = %10.6f min^-1' % get_value(Rx_m))
print('u(Rx)      = %10.6f min^-1' % get_stdunc(Rx_m))
print('Rs         = %10.6f min^-1' % get_value(Rs_m))
print('u(Rs)      = %10.6f min^-1' % get_stdunc(Rs_m))
print('r(Rx,Rs)   = %10.6f' % corr[0, 1])

R_m = Rx_m / Rs_m
print('R          = %10.6f' % get_value(R_m))
print('u(R)       = %10.6f' % get_stdunc(R_m))

# Calculation of final results
As = ufloat(0.1368, 0.0018) # Bq/g
ms = ufloat(5.0192, 0.0050) # g
mx = ufloat(5.0571, 0.0010) # g

Ax = As * ms / mx * R_m
print('Final result:')
print('Ax         = %10.6f Bq/g' % get_value(Ax))
print('u(Ax)      = %10.6f Bq/g' % get_stdunc(Ax))
print('u(Ax)/Ax   = %10.6f Bq/g' % (get_stdunc(Ax) / get_value(Ax)))

print('End')
