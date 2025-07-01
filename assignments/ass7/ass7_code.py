# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 14:49:29 2024

@author: polop
"""

import numpy as np
from scipy.integrate import solve_ivp

# Function f corresponding to the ODE system
def f(t, x, mu, dir):
    df = np.zeros(4 + 4*4)

    r1 = np.sqrt((x[0] - mu)**2 + x[1]**2)
    r2 = np.sqrt((x[0] - mu + 1)**2 + x[1]**2)

    # Equations for position and velocity
    df[0] = x[2]
    df[1] = x[3]
    df[2] = 2 * x[3] + x[0] - ((1 - mu) * (x[0] - mu) / r1**3) - mu * (x[0] - mu + 1) / r2**3
    df[3] = -2 * x[2] + x[1] * (1 - (1 - mu) / r1**3 - mu / r2**3)

    # Omegaxx, Omegayy, Omegaxy terms
    Omegaxx = (mu - 1) / ((mu - x[0])**2 + x[1]**2)**(3/2) - mu / ((x[0] - mu + 1)**2 + x[1]**2)**(3/2) \
              - (3 * (2 * mu - 2 * x[0])**2 * (mu - 1)) / (4 * ((mu - x[0])**2 + x[1]**2)**(5/2)) \
              + (3 * mu * (2 * x[0] - 2 * mu + 2)**2) / (4 * ((x[0] - mu + 1)**2 + x[1]**2)**(5/2)) + 1

    Omegayy = (mu - 1) / ((mu - x[0])**2 + x[1]**2)**(3/2) - mu / ((x[0] - mu + 1)**2 + x[1]**2)**(3/2) \
              + (3 * mu * x[1]**2) / ((x[0] - mu + 1)**2 + x[1]**2)**(5/2) \
              - (3 * x[1]**2 * (mu - 1)) / ((mu - x[0])**2 + x[1]**2)**(5/2) + 1

    Omegaxy = (3 * x[1] * (2 * mu - 2 * x[0]) * (mu - 1)) / (2 * ((mu - x[0])**2 + x[1]**2)**(5/2)) \
              + (3 * mu * x[1] * (2 * x[0] - 2 * mu + 2)) / (2 * ((x[0] - mu + 1)**2 + x[1]**2)**(5/2))

    # Variational equations for the monodromy matrix
    df[4:8] = x[12:16]
    df[8:12] = x[16:20]

    df[12] = Omegaxx * x[4] + Omegaxy * x[8] + 2 * x[16]
    df[13] = Omegaxx * x[5] + Omegaxy * x[9] + 2 * x[17]
    df[14] = Omegaxx * x[6] + Omegaxy * x[10] + 2 * x[18]
    df[15] = Omegaxx * x[7] + Omegaxy * x[11] + 2 * x[19]

    df[16] = Omegaxy * x[4] + Omegayy * x[8] - 2 * x[12]
    df[17] = Omegaxy * x[5] + Omegayy * x[9] - 2 * x[13]
    df[18] = Omegaxy * x[6] + Omegayy * x[10] - 2 * x[14]
    df[19] = Omegaxy * x[7] + Omegayy * x[11] - 2 * x[15]

    if dir == -1:
        df = -df

    return df

# Parameters
mu = 0.012
dir = -1
tol = 1e-13

# Initial conditions
x0 = [1.119748482176185, 0, 0, -0.2260036463376957]
x0 = np.concatenate((x0, np.reshape(np.eye(4), 16)))  # Monodromy matrix in flat form

# Time span
tf = 1 * 6.219656012174353

# Solve ODE system
sol = solve_ivp(f, [0, tf], x0, args=(mu, dir), atol=tol, rtol=tol)

# Compute determinants of the 4x4 matrix at each time step
dets = []
for i in range(len(sol.y[0])):
    A = np.array([sol.y[4:8, i], sol.y[8:12, i], sol.y[12:16, i], sol.y[16:20, i]])
    dets.append(np.linalg.det(A))

# Difference between max and min determinants
determinant_difference = np.max(dets) - np.min(dets)
print(determinant_difference)

import matplotlib.pyplot as plt

# Gráfico de Y[:,1] (posición X) en función del tiempo
plt.figure(figsize=(8, 5))
plt.plot(sol.y[0], sol.y[1], label="X (posición en el eje X)", color='blue')
plt.xlabel("X")
plt.ylabel("Y")
plt.title("phase portrait")
plt.legend()
plt.grid()
plt.show()

