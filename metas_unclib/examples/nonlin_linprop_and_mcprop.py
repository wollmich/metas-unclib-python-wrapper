# Example Non-Linear using LinProp and MCProp
# Michael Wollensack METAS - 22.02.2023

print('Example Non-Linear using LinProp and MCProp')
print('Begin')

from metas_unclib import *
use_linprop()

# Definition of the inputs using LinProp
x = ufloat(0.0, 1.0, desc='x')

# Compute the outputs using LinProp
y = x * x
z = 10 * y

outputs = [y, z]
cv = get_covariance(outputs)
print(outputs)
print(cv)

# Convert LinProp to MCProp
xMC = uspecial.linprop2mcprop(x)

# Compute the outputs using MCProp
yMC = xMC * xMC
zMC = 10 * yMC

outputsMC = [yMC, zMC]
cvMC = get_covariance(outputsMC)
print(outputsMC)
print(cvMC)

# Convert MCProp to LinProp
outputs2 = uspecial.mcprop2linprop(outputsMC, xMC, x)
cv2 = get_covariance(outputs2)
print(outputs2)
print(cv2)

# Test
test = cv2 - cvMC
print(test)

print('End')
