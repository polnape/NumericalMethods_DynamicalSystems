# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 21:45:29 2024

@author: polop
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
# Ejecución del código
mu = 0.1
tol = 1e-13
index = 1
stable = -1
s = 1e-6
iregion = 1
# compute Li given mu
# return functions to find epsilon given a mu value
def quintic1(xi, mu):
    """Función quintic1 para el cálculo de L1."""
    return ((mu * (1 - xi) ** 2) / (3 - 2 * mu - xi * (3 - mu - xi))) ** (1 / 3)

def quintic2(xi, mu):
    """Función quintic2 para el cálculo de L2."""
    return ((mu * (1 + xi) ** 2) / (3 - 2 * mu + xi * (3 - mu + xi))) ** (1 / 3)

def quintic3(xi, mu):
    """Función quintic3 para el cálculo de L3."""
    return (((1 - mu) * (1 + xi) ** 2) / (1 + 2 * mu + xi * (2 + mu + xi))) ** (1 / 3)

def Omega(x, mu, r1, r2):

    return 0.5 * (x[0]**2 + x[1]**2) + (1 - mu) / r1 + mu / r2 + 0.5 * mu * (1 - mu)

def Omegaxx(x, y, mu):
    """ Calcula Ω_xx dado x, y, y mu. """
    r1 = np.sqrt((x - mu)**2 + y**2)
    r2 = np.sqrt((x - mu + 1)**2 + y**2)
    term1 = 1 - (1 - mu) / r1**3
    term2 = (3 * (1 - mu) * (x - mu)**2) / r1**5
    term3 = -mu / r2**3
    term4 = (3 * mu * (x - mu + 1)**2) / r2**5
    return term1 + term2 + term3 + term4

def Omegaxy(x, y, mu):
    """ Calcula Ω_xy dado x, y, y mu. """
    r1 = np.sqrt((x - mu)**2 + y**2)
    r2 = np.sqrt((x - mu + 1)**2 + y**2)
    term1 = (3 * (1 - mu) * (x - mu) * y) / r1**5
    term2 = (3 * mu * y) / r2**5
    return y * (term1 + term2)

def Omegayy(x, y, mu):
    """ Calcula Ω_yy dado x, y, y mu. """
    r1 = np.sqrt((x - mu)**2 + y**2)
    r2 = np.sqrt((x - mu + 1)**2 + y**2)
    term1 = 1 - (1 - mu) / r1**3 - mu / r2**3
    term2 = y**2 * ((3 * (1 - mu)) / r1**5 + (3 * mu) / r2**5)
    return term1 + term2


def itera(mu, x0, f, tol):
    xi = x0
    xik = f(xi, mu)
    while(abs(xi-xik)>tol):
        xi = xik
        xik = f(xik, mu)
   
    xi = xik
    return xi

x0 = (mu/(3*(1-mu)))**(1/3)
# Computar L1, L2, o L3 según el valor de index

if index == 1:
    xi = itera(mu, (mu / (3 * (1 - mu)))**(1 / 3), quintic1, tol)
    Lix = mu - 1 + xi
    C = 2*Omega([Lix, 0], mu, 1-xi, xi)
elif index == 2:
    xi = itera(mu, (mu / (3 * (1 - mu)))**(1 / 3), quintic2, tol)
    Lix = mu - 1 - xi
    C = 2 * Omega([Lix,0], mu, np.sqrt((Lix - mu)**2 + 0**2), abs(Lix - mu + 1))
elif index == 3:
    xi = itera(mu, 1 - 7 / 12 * mu, quintic3, tol)
    Lix = mu + xi
    C = 2*Omega([Lix,0], mu, xi, 1+xi)

print(f"Valor de L{index}x = {Lix} and the Jacobi Constant is C = {C}")





# Cálculo de la matriz Jacobiana
def jacobian_matrix(Lix, mu):
    A = np.zeros((4, 4))
    A[0, 2] = 1
    A[1, 3] = 1
    A[2, 3] = 2
    A[3, 2] = -2
    A[2, 0] = Omegaxx(Lix, 0, mu)
    A[2, 1] = Omegaxy(Lix, 0, mu)
    A[3, 0] = Omegaxy(Lix, 0, mu)
    A[3, 1] = Omegayy(Lix, 0, mu)
    return A


import numpy as np
from scipy.linalg import eig
def compute_eigenvectors(A, stable):
    """
    Calcula y selecciona el eigenvector real basado en el signo del eigenvalor.
    
    Parámetros:
    A - Matriz cuadrada (numpy array).
    stable - 1 para manifold estable (eigenvalor negativo),
             -1 para manifold inestable (eigenvalor positivo).
    
    Retorna:
    Eigenvector real correspondiente al criterio seleccionado como un array de tipo float.
    """
    # Calcular todos los eigenvalores y eigenvectores usando scipy.linalg.eig
    eigvals, eigvecs = eig(A)
    
    # Inicializar la variable para almacenar el índice del eigenvector correcto
    c = None
    
    # Iterar sobre todos los eigenvalores para seleccionar el correcto según el criterio
    for i in range(len(eigvals)):
        if stable == 1 and np.real(eigvals[i]) < 0 and np.isclose(np.imag(eigvals[i]), 0, atol=1e-9):
            c = i
        elif stable == -1 and np.real(eigvals[i]) > 0 and np.isclose(np.imag(eigvals[i]), 0, atol=1e-9):
            c = i
    
    # Verificar si se encontró un eigenvalor adecuado
    if c is not None:
        # Extraer el eigenvector correspondiente y convertir a real si la parte imaginaria es insignificante
        eigenvector_real = np.real(eigvecs[:, c])
        return eigenvector_real.astype(float)
    else:
        raise ValueError("No se encontró un eigenvalor adecuado para el criterio dado.")


# Simulación de las variedades inestables/estables
def compute_manifold(Li, v, s, mu, iregion):
    # Ajustar el signo del vector según iregion
    if iregion > 0 and v[1] < 0:
        v = -v
    elif iregion < 0 and v[1] > 0:
        v = -v

    # Condición inicial ajustada
    x0 = Li[0] + s * v[0]
    y0 = Li[1] + s * v[1]
    vx0 = s * v[2]
    vy0 = s * v[3]
    
    return [x0, y0, vx0, vy0]



# Cálculo de la matriz Jacobiana y los eigenvectores
A = jacobian_matrix(Lix, mu)
try:
    eigenvector_inestable = compute_eigenvectors(A, stable=-1)
    print("Eigenvector inestable:", eigenvector_inestable)
except ValueError as e:
    print(e)

# Calcular eigenvector para el manifold estable
try:
    eigenvector_estable = compute_eigenvectors(A, stable=1)
    print("Eigenvector estable:", eigenvector_estable)
except ValueError as e:
    print(e)
#%%
# In case if we need it
def newton_raphson(f, Df, x0, tol, a):
    xk = np.array(x0, dtype=float)  # Inicialización del punto
    tolk = 10.0  # Inicialización de la tolerancia

    while tolk > tol:
        # Calcular la corrección dxk usando la derivada de f y el valor actual xk
        # término -f(xk)/(f'(xk))
        dxk = np.linalg.solve(Df(xk, a), -f(xk, a))
        # nuevo xk, me acerco a la raíz
        xk = xk + dxk

        # Calcular la norma del ajuste en L2, estamos usando una tolerancia horizontal
        tolk = np.linalg.norm(dxk, 2)

        # 
        # Ponemos el nuevo punto entre 0 y 2pi
        xk[0] = np.mod(xk[0], 2 * np.pi)
        xk[1] = np.mod(xk[1], 2 * np.pi)

    return xk

def g(pos):
    return pos[1]
def grad(g):
    dg = np.array([0,1])
    return dg

def Poincar_Map(initial_point, ncrossings, idir, function_to_evaluate):
    """
    Calcula la trayectoria completa hasta que se crucen la sección Poincaré (y = 0) dos veces.
    La integración se detiene una vez que se alcanzan los cruces especificados.
    """
    crossings = []
    results_x = []
    results_y = []

    # Definir la dirección del tiempo
    t_periodic = [0, idir * 100]
    
    # Realizar la integración hasta alcanzar el número deseado de cruces
    sol = solve_ivp(function_to_evaluate, t_periodic, initial_point, t_eval=np.linspace(0, idir * 100, 10000), rtol=1e-9, atol=1e-12)
    
    # Iterar sobre los puntos integrados para detectar los cruces y almacenar la trayectoria
    for i in range(1, len(sol.t)):
        y1 = sol.y[1, i-1]
        y2 = sol.y[1, i]

        # Almacenar la trayectoria completa
        results_x.append(sol.y[0, i])
        results_y.append(sol.y[1, i])

        # Verificar si hay un cruce (cambio de signo en y)
        if y1 * y2 < 0:
            # Interpolación lineal para encontrar el cruce exacto
            x_before, y_before = sol.y[0, i-1], sol.y[1, i-1]
            x_after, y_after = sol.y[0, i], sol.y[1, i]
            t_before, t_after = sol.t[i-1], sol.t[i]

            # Interpolación para el cruce
            cross_t = t_before + (t_after - t_before) * (-y_before) / (y_after - y_before)
            cross_x = x_before + (x_after - x_before) * (-y_before) / (y_after - y_before)
            
            # Agregar el cruce encontrado
            crossings.append([cross_x, 0])
            
            # Detener la integración si alcanzamos el número deseado de cruces
            if len(crossings) >= ncrossings:
                break

    return results_x, results_y, crossings






# sitema de equaciones
def f(t, x, mu, dir):

    # Inicializar el vector de derivadas
    df = np.zeros(4)
    
    # Calcular las distancias r1 y r2
    r1 = np.sqrt((x[0] - mu)**2 + x[1]**2)
    r2 = np.sqrt((x[0] - mu + 1)**2 + x[1]**2)
    
    # Definir el sistema de ecuaciones diferenciales
    df[0] = x[2]  # dx/dt = v1
    df[1] = x[3]  # dy/dt = v2
    df[2] = 2 * x[3] + x[0] - ((1 - mu) * (x[0] - mu) / r1**3) - (mu * (x[0] - mu + 1) / r2**3)
    df[3] = -2 * x[2] + x[1] * (1 - (1 - mu) / r1**3 - mu / r2**3)
    
    # Ajustar la dirección de integración si dir == -1
    if dir == -1:
        df = -df
    
    return df




# Definir condiciones iniciales para ambas ramas
x1 = [Lix, 0, 0, 0] + s * eigenvector_inestable
x2 = [Lix, 0, 0, 0] - s * eigenvector_inestable





ncrossings = 2  # Número de cruces deseado en la sección Poincaré

# EN LA FUNCIÓN DE POINCARÉ YA ME ENCARGO DE GRAFICAR TODO, HASTA LLEGAR AL NÚMERO DE CROSSINGS COMPLETO
# Calcular los cruces y las trayectorias para la rama positiva
results_x1, results_y1, crossings_positive = Poincar_Map(x1, ncrossings, idir=1, function_to_evaluate=lambda t, x: f(t, x, mu, 1))

# Calcular los cruces y las trayectorias para la rama negativa
results_x2, results_y2, crossings_negative = Poincar_Map(x2, ncrossings, idir=1, function_to_evaluate=lambda t, x: f(t, x, mu, 1))


# Graficar las trayectorias completas para ambas ramas
plt.figure()
plt.plot(results_x1, results_y1, c='orange', label="Unstable Manifold (positive branch)")
plt.plot(results_x2, results_y2, c='blue', label="Unstable Manifold (negative branch)")

# Graficar los cruces detectados
plt.scatter([p[0] for p in crossings_positive], [p[1] for p in crossings_positive], c='red')
plt.scatter([p[0] for p in crossings_negative], [p[1] for p in crossings_negative], c='green')

plt.xlabel("x")
plt.ylabel("y")
plt.legend(loc = "upper left")
plt.grid()
plt.title("Manifolds Inestables hasta el segundo cruce")
plt.show()

#%%
x1_estable = [Lix, 0, 0, 0] + s * eigenvector_estable
x2_estable = [Lix, 0, 0, 0] - s * eigenvector_estable
# Definir condiciones iniciales para ambas ramas
x1 = [Lix, 0, 0, 0] + s * eigenvector_inestable
x2 = [Lix, 0, 0, 0] - s * eigenvector_inestable


ncrossings = 2  # Número de cruces deseado en la sección Poincaré

# EN LA FUNCIÓN DE POINCARÉ YA ME ENCARGO DE GRAFICAR TODO, HASTA LLEGAR AL NÚMERO DE CROSSINGS COMPLETO
# Calcular los cruces y las trayectorias para la rama positiva
results_x1, results_y1, crossings_positive = Poincar_Map(x1, ncrossings, idir=1, function_to_evaluate=lambda t, x: f(t, x, mu, 1))

# Calcular los cruces y las trayectorias para la rama negativa
results_x2, results_y2, crossings_negative = Poincar_Map(x2, ncrossings, idir=1, function_to_evaluate=lambda t, x: f(t, x, mu, 1))
results_x1_es, results_y1_es, crossings_positive = Poincar_Map(x1_estable, ncrossings, idir=1, function_to_evaluate=lambda t, x: f(t, x, mu, -1))

# Calcular los cruces y las trayectorias para la rama negativa
results_x2_es, results_y2_es, crossings_negative = Poincar_Map(x2_estable, ncrossings, idir=1, function_to_evaluate=lambda t, x: f(t, x, mu, -1))


# Graficar las trayectorias completas para ambas ramas
plt.figure()
plt.plot(results_x1, results_y1, c='orange', label="Unstable Manifold (positive branch)")
plt.plot(results_x2, results_y2, c='blue', label="Unstable Manifold (negative branch)")
plt.plot(results_x1_es, results_y1_es, c='orange', label="Stable Manifold (positive branch)")
plt.plot(results_x2_es, results_y2_es, c='blue', label="Stable Manifold (negative branch)")

# Graficar los cruces detectados
plt.scatter([p[0] for p in crossings_positive], [p[1] for p in crossings_positive], c='red')
plt.scatter([p[0] for p in crossings_negative], [p[1] for p in crossings_negative], c='green')

plt.xlabel("x")
plt.ylabel("y")
plt.legend(loc = "upper left")
plt.grid()
plt.title("Manifolds Inestables hasta el segundo cruce")
plt.show()