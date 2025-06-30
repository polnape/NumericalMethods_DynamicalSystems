# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 17:22:08 2024

@author: polop
"""

#%%
# PART 1: LOGISTIC MAP
# valores de las iteraciones (cada iteración tiene dos puntos x0 y xf)
    
import numpy as np
import matplotlib.pyplot as plt

# Parámetros
a = 3.2
# a = 2.0
x0 = 0.4
# x0 = 0.2

n = 100  # Número de iteraciones

# Definir el mapa logístico
def logistic_map(x, a):
    return a * x * (1 - x)

# Parábola
x_parabola = np.linspace(0, 1, 500)
parabola = logistic_map(x_parabola, a)

# Inicializar el array para las iteraciones
iterations = np.zeros((2 * n, 2))
iterations[0, 0] = x0
iterations[0, 1] = 0

# Iterar para llenar el array
i = 0
x = x0
y = 0
while i < 2 * n - 1:
    iterations[i, 0] = x
    iterations[i, 1] = y
    y = logistic_map(x, a)
    i += 1
    iterations[i, 0] = x
    iterations[i, 1] = y
    x = y
    i += 1
 


# Graficar
plt.figure(figsize=(10, 6))

# Graficar la parábola
plt.plot(x_parabola, parabola, label=f'y = {a} * x * (1 - x)', linewidth=1.5)

# Graficar la recta x = y
plt.plot(x_parabola, x_parabola, label='y = x', linestyle='--', linewidth=1.5, color = "green")

# Graficar los puntos de iteraciones
plt.plot(iterations[:, 0], iterations[:, 1], label='Iteraciones', color='red', marker='o', markersize=3)

# Añadir etiquetas y título
plt.xlabel('x')
plt.ylabel('y')
# plt.title('Parábola y = a * x * (1 - x) con Iteraciones')
plt.legend()
plt.ylim
plt.show()

#%%
   
# PART 2: STANDARD MAP

import numpy as np
import matplotlib.pyplot as plt


a = -0.7
NP = 100  # Número de condiciones iniciales
iteraciones = 500  # Número de iteraciones

xf = np.zeros((iteraciones, 2))

IC = np.linspace(-np.pi, np.pi, NP)  

def f(x, a):
    # valor de la función F(x,y)
    F_xy = np.zeros(2)
    # valor x y lo envolvemos en el rango de 0 a 2pi
    F_xy[0] = np.mod(x[0] + a * np.sin(x[0] + x[1]), 2 * np.pi)
    F_xy[1] = np.mod(x[0] + x[1], 2 * np.pi)

    # pasamos los puntos a [-pi, pi]
    if F_xy[0] > np.pi and F_xy[1] > np.pi:
        F_xy[0] = -np.pi + np.mod(F_xy[0], np.pi)
        F_xy[1] = -np.pi + np.mod(F_xy[1], np.pi)

    if F_xy[0] > np.pi and F_xy[1] < np.pi:
        F_xy[0] = -np.pi + np.mod(F_xy[0], np.pi)

    if F_xy[0] < np.pi and F_xy[1] > np.pi:
        F_xy[1] = -np.pi + np.mod(F_xy[1], np.pi)

    return F_xy


# Crear un mapa de colores
cmap = plt.cm.get_cmap('hsv', NP)  # Utilizamos 'hsv' para obtener una amplia gama de colores
# Configurar el gráfico
plt.figure(figsize=(8, 8))
plt.xlim([-np.pi, np.pi])
plt.ylim([-np.pi, np.pi])
plt.title("Standard 2D map for a = -0.7")

# Bucle para todas las initial conditions que tenemos
for i in range(NP):
    x0 = np.array([IC[i], 0])
    xf[0, :] = f(x0, a)
    for j in range(1, iteraciones):
        xf[j, :] = f(xf[j-1, :], a)
    
    # Asignar un color diferente a cada condición inicial usando cmap
    color = cmap(i)
    plt.scatter(xf[:, 0], xf[:, 1], s=10, color=color)

# # Marcar el 2-periodic point con una cruz negra
# plt.scatter(2.7967814286635737, 1.7432019392580074, color='black', marker='x', s=100, label="2-periodic point")

# # Marcar el 3-periodic point con una cruz marron
# plt.scatter(1.7380081367622715, 2.2725885852086543, color='brown', marker='x', s=100, label="3-periodic point")
# plt.legend()
plt.show()

#%%
# Find exact initial condition of a 2-periodic point
# Después de iterar dos veces obtenemos el mismo punto
# f(f(x0,a), a) = x0

#newton_raphson method, para encontrar raices
def newton_raphson(f, Df, x0, tolerance, a):
    xk = np.array(x0, dtype=float)  # Inicialización del punto
    current_tol = 10.0  # Inicialización de la tolerancia

    while current_tol > tolerance:
        # Calcular la corrección dxk usando la derivada de f y el valor actual xk
        # término -f(xk)/(f'(xk))
        dxk = np.linalg.solve(Df(xk, a), -f(xk, a))
        # nuevo xk, me acerco a la raíz
        xk = xk + dxk

        # Calcular la norma del ajuste en L2, estamos usando una tolerancia horizontal
        current_tol = np.linalg.norm(dxk, 2)

        # 
        # Ponemos el nuevo punto entre 0 y 2pi
        xk[0] = np.mod(xk[0], 2 * np.pi)
        xk[1] = np.mod(xk[1], 2 * np.pi)

    return xk




# Definir el mapa estándar (función f)
def f(x, a):
    f_xy = np.zeros(2)
    f_xy[0] = np.mod(x[0] + a * np.sin(x[0] + x[1]), 2 * np.pi)
    f_xy[1] = np.mod(x[0] + x[1], 2 * np.pi)

    # Ajustar los puntos de [0, 2*pi] a [-pi, pi]
    if f_xy[0] > np.pi:
        f_xy[0] -= 2 * np.pi
    if f_xy[1] > np.pi:
        f_xy[1] -= 2 * np.pi

    return f_xy

# Matriz Jacobiana, de las respectivas derivadas
def Df(x, a):
    Df_matrix = np.zeros((2, 2))
    Df_matrix[0, 0] = 1 + a * np.cos(x[0] + x[1])
    Df_matrix[0, 1] = a * np.cos(x[0] + x[1])
    Df_matrix[1, 0] = 1
    Df_matrix[1, 1] = 1
    return Df_matrix


# Defino la función a minimizar minf = f^k(x,a) - x para un k periodic point
# Función a minimizar
def minf(k, x, a): 
    result = f(x, a)
    for _ in range(1, k):
        result = f(result, a)
    return result - x

# Jacobiano de la función a minimizar
def Derminf(k, x, a):
    jacobiano = Df(x, a)  # Inicia el Jacobiano con la primera derivada Df(x, a)
    x_copy = np.copy(x)    # Copia de x para no modificar el valor original
    for _ in range(1, k):
        x_copy = f(x_copy, a)  # Aplica f para obtener el siguiente valor de x
        jacobiano = np.dot(Df(x_copy, a), jacobiano)  # Calcula el Jacobiano iterativamente
    return jacobiano - np.eye(2)  # Resta la matriz identidad al final

#%%
# Parámetros iniciales
a = -0.7
tolerancia = 1e-15 
x0 = np.array([0.5, 0.5])  
# epsilon le he dado un valor muy pequeño porque si no se me desviava mucho de lo esperado
epsilon = np.array([0.02, 0.02]) 
# k periodic point
k = int(input("Introduce el valor entero de k: "))


# Usar lambdas para no ejecutar las funciones inmediatamente
per_point = newton_raphson(lambda xk, a: minf(k, xk, a), lambda xk, a: Derminf(k, xk, a), x0, tolerancia, a)
orbit = np.zeros((k,2))

# Mostrar el punto periódico encontrado
print(f'The {k} periodic point is located at ({per_point[0]}, {per_point[1]})')


# Plot the associated orbit. Plot an invariant curve around it

# Construir la órbita del movimiento en el punto periódico después de las iteraciones
orbit[0, :] = per_point[0]
# Iteramos las k veces para volver al punto
for i in range(1, k):
    # valores de los puntos de la órbita
    orbit[i, :] = f(orbit[i - 1, :], a)

# nos desplazamos un epsilon del k-periodic point para ver como es la curva invariante y estudiar el k-periodic point
x = per_point[0] + epsilon
# n (iteraciones) puntos, bidimensional
n = 1000
curve = np.zeros((n, 2))
# empezamos en el k periodic point
curve[0, :] = x
for i in range(1, n):
    curve[i, :] = f(curve[i - 1, :], a)

# Graficar la órbita
plt.figure(figsize=(8, 8))
plt.scatter(curve[:,0], curve[:,1], s=20, color='blue', label='Invariant curve')
# la órbita tendrá k puntos
plt.scatter(orbit[:, 0], orbit[:, 1], s=30, color='orange', label='Periodic orbit')
plt.xlim([-np.pi, np.pi])
plt.ylim([-np.pi, np.pi])
plt.title(f"2D standard map orbit for k={k}, a={a}")
plt.xlabel("x")
plt.ylabel("y")
plt.scatter(per_point[0], per_point[1], color='red', marker='x', s=100, label=f'{k} periodic point')
plt.legend()

plt.show()

#%%
 # take a=-0.1,-0.3,-0.5,.....,-2.1
 # (10 different values of a), and for each a, obtain the
 # output data. Plot, as a film, the evolution of the dynamics
 # varying a, that is, plot 10 different plots.


a_list = np.arange(-0.1, -2.2, -0.2, dtype = float)
for el in a_list:
    # Crear un mapa de colores
    cmap = plt.cm.get_cmap('hsv', NP)  # Utilizamos 'hsv' para obtener una amplia gama de colores
    # Configurar el gráfico
    plt.figure(figsize=(8, 8))
    plt.xlim([-np.pi, np.pi])
    plt.ylim([-np.pi, np.pi])
    plt.title(f"Standard 2D map for a = {el:.1f}")

    # Bucle para todas las initial conditions que tenemos
    for i in range(NP):
        x0 = np.array([IC[i], 0])
        xf[0, :] = f(x0, el)
        for j in range(1, iteraciones):
            xf[j, :] = f(xf[j-1, :], el)
        
        # Asignar un color diferente a cada condición inicial usando cmap
        color = cmap(i)
        plt.scatter(xf[:, 0], xf[:, 1], s=10, color=color)


    plt.show()
#%%
# 2D MAP with 2 and 3 periodic points


# Crear un mapa de colores
cmap = plt.cm.get_cmap('hsv', NP)  # Utilizamos 'hsv' para obtener una amplia gama de colores
# Configurar el gráfico
plt.figure(figsize=(8, 8))
plt.xlim([-np.pi, np.pi])
plt.ylim([-np.pi, np.pi])
plt.title("Standard 2D map for a = -0.7")

# Bucle para todas las initial conditions que tenemos
for i in range(NP):
    x0 = np.array([IC[i], 0])
    xf[0, :] = f(x0, a)
    for j in range(1, iteraciones):
        xf[j, :] = f(xf[j-1, :], a)
    
    # Asignar un color diferente a cada condición inicial usando cmap
    color = cmap(i)
    plt.scatter(xf[:, 0], xf[:, 1], s=10, color="black")

# Marcar el 2-periodic point con una cruz negra
plt.scatter(2.7967814286635737, 1.7432019392580074, color='green', marker='x', s=100, label="2-periodic point")

# Marcar el 3-periodic point con una cruz marron
plt.scatter(1.7380081367622715, 2.2725885852086543, color='red', marker='x', s=100, label="3-periodic point")



plt.legend()
plt.show()

#%%

#  Plot a vertical invariant curve. Plot another one around the origin.

# Parámetros
n = 1000  # Número de puntos en la curva
epsilon = np.array([0.01, 0.0])  # Pequeño desplazamiento para la curva vertical

# Curva invariante vertical
vertical_invariant_curve = np.zeros((n, 2))
# La curva invariante vertical empezamos en un punto que queremos estudiar y lo vamos perturbando con x constante, solo movimento la y
vertical_invariant_curve[0, :] = [2, 0 + epsilon[1]]
for i in range(1, n):
    vertical_invariant_curve[i, :] = f(vertical_invariant_curve[i - 1, :], a)

# Curva invariante en el origen
origin_curve = np.zeros((n, 2))
# La curva invariante en el origen empiezo en un punto x arbritario y = 0, y lo voy perturbando en x
origin_curve[0, :] = [1 + epsilon[0], 0]
for i in range(1, n):
    origin_curve[i, :] = f(origin_curve[i - 1, :], a)

# Graficar las curvas invariantes junto con las órbitas
plt.figure(figsize=(8, 8))

plt.scatter(vertical_invariant_curve[:, 0], vertical_invariant_curve[:, 1], s=10, color='green', label='Curva invariante vertical')
plt.scatter(origin_curve[:, 0], origin_curve[:, 1], s=10, color='purple', label='Curva invariante en el origen')

plt.xlim([-np.pi, np.pi])
plt.ylim([-np.pi, np.pi])
plt.title(f"Curvas invariantes y órbitas para k={k}, a={a}")
plt.xlabel("x")
plt.ylabel("y")

plt.legend()
plt.show()
