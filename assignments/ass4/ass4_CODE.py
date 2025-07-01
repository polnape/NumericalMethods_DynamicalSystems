# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 15:55:44 2024

@author: polop
"""


import numpy as np
import matplotlib.pyplot as plt

def itera(mu, x0, f, tol):
    xi = x0
    xik = f(xi, mu)
    while abs(xi - xik) > tol:
        xi = xik
        xik = f(xik, mu)
    return xik

def quintic1(xi, mu):
    return ((mu * (1 - xi) ** 2) / (3 - 2 * mu - xi * (3 - mu - xi))) ** (1 / 3)

def quintic2(xi, mu):
    return ((mu * (1 + xi) ** 2) / (3 - 2 * mu + xi * (3 - mu + xi))) ** (1 / 3)

def quintic3(xi, mu):
    return (((1 - mu) * (1 + xi) ** 2) / (1 + 2 * mu + xi * (2 + mu + xi))) ** (1 / 3)

def Omega(x, mu, r1, r2):
    return 0.5 * (x[0] ** 2 + x[1] ** 2) + (1 - mu) / r1 + mu / r2 + 0.5 * (mu * (1 - mu))

def Omegaxx(x, mu):
    term1 = (mu - 1) / ((mu - x[0]) ** 2 + x[1] ** 2) ** (3 / 2)
    term2 = -mu / ((x[0] - mu + 1) ** 2 + x[1] ** 2) ** (3 / 2)
    term3 = -(3 * (2 * mu - 2 * x[0]) ** 2 * (mu - 1)) / (4 * ((mu - x[0]) ** 2 + x[1] ** 2) ** (5 / 2))
    term4 = (3 * mu * (2 * x[0] - 2 * mu + 2) ** 2) / (4 * ((x[0] - mu + 1) ** 2 + x[1] ** 2) ** (5 / 2))
    return term1 + term2 + term3 + term4 + 1

def Omegayy(x, mu):
    term1 = (mu - 1) / ((mu - x[0]) ** 2 + x[1] ** 2) ** (3 / 2)
    term2 = -mu / ((x[0] - mu + 1) ** 2 + x[1] ** 2) ** (3 / 2)
    term3 = (3 * mu * x[1] ** 2) / ((x[0] - mu + 1) ** 2 + x[1] ** 2) ** (5 / 2)
    term4 = -(3 * x[1] ** 2 * (mu - 1)) / ((mu - x[0]) ** 2 + x[1] ** 2) ** (5 / 2)
    return term1 + term2 + term3 + term4 + 1

def Omegaxy(x, mu):
    term1 = (3 * x[1] * (2 * mu - 2 * x[0]) * (mu - 1)) / (2 * ((mu - x[0]) ** 2 + x[1] ** 2) ** (5 / 2))
    term2 = (3 * mu * x[1] * (2 * x[0] - 2 * mu + 2)) / (2 * ((x[0] - mu + 1) ** 2 + x[1] ** 2) ** (5 / 2))
    return term1 + term2

# Parameters
tol = 1e-14
mu = 0.008
index = 3
n = 1000

# Computation of Li
if index == 1:
    xi = itera(mu, (mu / (3 * (1 - mu))) ** (1 / 3), quintic1, tol)
    L1x = mu - 1 + xi
    print(L1x)
    print(2 * Omega([L1x, 0], mu, 1 - xi, xi))
elif index == 2:
    xi = itera(mu, (mu / (3 * (1 - mu))) ** (1 / 3), quintic2, tol)
    L2x = mu - 1 - xi
    print(L2x)
    print(2 * Omega([L2x, 0], mu, ((L2x - mu) ** 2 + 0 ** 2) ** (1 / 2), abs(L2x - mu + 1)))
elif index == 3:
    xi = itera(mu, 1 - 7 / 12 * mu, quintic3, tol)
    L3x = mu + xi
    print(L3x)
    print(2 * Omega([L3x, 0], mu, xi, 1 + xi))

# Plot (mu, x), (mu, C) and (mu, char. exp.)
mu_vec = np.linspace(0, 0.5, n)
xL1, xL2, xL3 = np.zeros(n), np.zeros(n), np.zeros(n)
CL1, CL2, CL3 = np.zeros(n), np.zeros(n), np.zeros(n)
ReL1, ImL1 = np.zeros((n, 4)), np.zeros((n, 4))

for i in range(n):
    xL1[i] = mu_vec[i] - 1 + itera(mu_vec[i], (mu_vec[i] / (3 * (1 - mu_vec[i]))) ** (1 / 3), quintic1, tol)
    xL2[i] = mu_vec[i] - 1 - itera(mu_vec[i], (mu_vec[i] / (3 * (1 - mu_vec[i]))) ** (1 / 3), quintic2, tol)
    xL3[i] = mu_vec[i] + itera(mu_vec[i], 1 - (7 / 12) * mu_vec[i], quintic3, tol)

    CL1[i] = 2 * Omega([xL1[i], 0], mu_vec[i], abs(xL1[i] - mu_vec[i]), abs(xL1[i] - mu_vec[i] + 1))
    CL2[i] = 2 * Omega([xL2[i], 0], mu_vec[i], abs(xL2[i] - mu_vec[i]), abs(xL2[i] - mu_vec[i] + 1))
    CL3[i] = 2 * Omega([xL3[i], 0], mu_vec[i], abs(xL3[i] - mu_vec[i]), abs(xL3[i] - mu_vec[i] + 1))

# Plot graphs
plt.figure(1)
plt.plot(mu_vec, xL1, label='xL1')
plt.plot(mu_vec, xL2, label='xL2')
plt.plot(mu_vec, xL3, label='xL3')
plt.legend()
plt.show()

plt.figure(2)
plt.plot(mu_vec, CL1, label='CL1')
plt.plot(mu_vec, CL2, label='CL2')
plt.plot(mu_vec, CL3, label='CL3')
plt.legend()
plt.show()

# P(xi, mu) for L1 to check that P has a unique zero.
xi_01 = np.linspace(0, 1, n)
P = np.zeros(n)
for i in range(n):
    P[i] = quintic1(xi_01[i], mu)

plt.figure()
plt.plot(xi_01, P)
plt.title('P(xi, mu) for L1')
plt.show()
 
