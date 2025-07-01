# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 14:51:45 2024

@author: polop
"""


#%%

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Definimos la función del sistema de ecuaciones diferenciales
def f(t, x, mu, direccion):
    r1 = np.sqrt((x[0] - mu)**2 + x[1]**2)
    r2 = np.sqrt((x[0] - mu + 1)**2 + x[1]**2)

    df = np.zeros(4)
    df[0] = x[2]
    df[1] = x[3]
    df[2] = 2 * x[3] + x[0] - ((1 - mu) * (x[0] - mu) / r1**3) - mu * (x[0] - mu + 1) / r2**3
    df[3] = -2 * x[2] + x[1] * (1 - (1 - mu) / r1**3 - mu / r2**3)

    if direccion == -1:
        df = -df
    # retorno x1, x2, x3, x4
    return df

# Omega definida en clases de teoría
def Omega(x, mu):
    r1 = np.sqrt((x[0] - mu)**2 + x[1]**2)
    r2 = np.sqrt((x[0] - mu + 1)**2 + x[1]**2)
    return (1/2) * (x[0]**2 + x[1]**2) + (1 - mu) / r1 + (mu / r2) + (1/2) * mu * (1 - mu)

# Parámetros del problema
tol = 1e-14
mu = 0.012
direccion = 1  # 1 para adelante en el tiempo, -1 para atrás
x0 = [1.119748482176185, 0, 0, -0.2260036463376957]
tf = 6.219656012174353

# Usamos solve_ivp para resolver el sistema de ecuaciones diferenciales
sol = solve_ivp(lambda t, x: f(t, x, mu, direccion), [0, tf], x0, rtol=tol, atol=tol, dense_output=True)

# Calculamos la constante de Jacobi en cada paso
JacobiCtt = np.zeros(len(sol.t))
for i in range(len(sol.t)):
    JacobiCtt[i] = 2 * Omega([sol.y[0, i], sol.y[1, i]], mu) - (sol.y[2, i]**2 + sol.y[3, i]**2)

# Graficamos la trayectoria
plt.plot(sol.y[0], sol.y[1])
plt.xlabel('x')
plt.ylabel('y')
plt.title('Órbita')
plt.grid(True)
plt.show()

# Graficamos la constante de Jacobi a lo largo del tiempo
plt.plot(sol.t, JacobiCtt)
plt.xlabel('Tiempo')
plt.ylabel('Constante de Jacobi')
plt.title('Conservación de la Constante de Jacobi')
plt.grid(True)
plt.show()
