# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 15:40:55 2024

@author: polop
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Parameters
tmax = [0, 2 * np.pi]
timestepctt = np.linspace(0, 100, 1000)

x0 = [1.0, 0.0]
tol = 1e-13
idir = -1

# ODE options for absolute and relative tolerance
opts = {'rtol': tol, 'atol': tol}

# Define the Hamiltonian function
def H(x):
    return 0.5 * (x[0]**2 + x[1]**2)

# Define the ODE system
def f(x, t, idir):
    dx = np.zeros(2)
    dx[0] = x[1]
    dx[1] = -x[0]
    if idir == -1:
        dx = -dx
    return dx

# Solve the ODE system using odeint
sol = odeint(f, x0, timestepctt, args=(idir,), **opts)

# Scatter plot of the solution on the plane
plt.scatter(sol[:, 0], sol[:, 1])
plt.xlabel('x')
plt.ylabel('y')
plt.title('Scatter plot of x(t) vs y(t)')
plt.show()
#%%
# Check if Hamiltonian is constant over time
max_H = -1e13
min_H = 1e13
for i in range(len(timestepctt)):
    eval_H = H(sol[i])
    if eval_H < min_H:
        min_H = eval_H
    if eval_H > max_H:
        max_H = eval_H
print("Hamiltonian Difference:", max_H - min_H)

# Variational equations
x0_VE = [1, 1, 1, 0, 0, 1]

# Define the variational equations
def VE(x, t, idir):
    dxVE = np.zeros(6)
    dxVE[0] = x[1]
    dxVE[1] = -x[0]
    dxVE[2] = x[4]
    dxVE[3] = x[5]
    dxVE[4] = -x[2]
    dxVE[5] = -x[3]
    if idir == -1:
        dxVE = -dxVE
    return dxVE

# Solve the variational equation
sol_VE = odeint(VE, x0_VE, timestepctt, args=(idir,), **opts)

# Check if the determinant is always 1
det = np.zeros(len(timestepctt))
for i in range(len(timestepctt)):
    det[i] = sol_VE[i, 2] * sol_VE[i, 5] - sol_VE[i, 3] * sol_VE[i, 4]

# Plot determinant over time to check if it's constant (equal to 1)
plt.plot(timestepctt, det)
plt.xlabel('t')
plt.ylabel('determinant')
plt.title('Determinant over time')

plt.show()
