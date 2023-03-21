# Example Right Triangle using LinProp and MCProp
# Michael Wollensack METAS - 22.02.2023

print('Example Right Triangle using LinProp and MCProp')
print('Begin')

from metas_unclib import *
use_linprop()

# Definition of the inputs using LinProp
a = ufloat(3.0, 0.3, desc='a')
b = ufloat(4.0, 0.4, desc='b')
inputs = [a, b]

# Compute the outputs using LinProp
c = umath.sqrt(a * a + b * b)
p = a + b + c
A = a * b / 2

outputs = [c, p, A]
cv = get_covariance(outputs)
print(outputs)
print(cv)

# Convert LinProp to MCProp
inputsMC = uspecial.linprop2mcprop(inputs)

# Compute the outputs using MCProp
aMC = inputsMC[0]
bMC = inputsMC[1]
cMC = umath.sqrt(aMC * aMC + bMC * bMC)
pMC = aMC + bMC + cMC
AMC = aMC * bMC / 2

outputsMC = [cMC, pMC, AMC]
cvMC = get_covariance(outputsMC)
print(outputsMC)
print(cvMC)

# Convert MCProp to LinProp
outputs2 = uspecial.mcprop2linprop(outputsMC, inputsMC, inputs)
cv2 = get_covariance(outputs2)
print(outputs2)
print(cv2)

# Test
test = cv2 - cvMC
print(test)

print('End')
