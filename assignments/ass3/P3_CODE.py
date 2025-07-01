# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 15:11:21 2024

@author: polop
"""
# PART A)



#%%

#%%

# PART A)

from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


#POINCARÉ SECTIONA AND POINCARÉ MAP

# PARTE A
# Poincaré section for y = 0

# The harmonic oscillator function
def oscilador_harmonico(t, pos):
    dx, dy = pos[1], -pos[0]  # Definir las derivadas
    return np.array([dx, dy])

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

# As we are working on y = 0, then x' = 0, the corresponding g(x) is defined by:
def g(pos):
    return pos[1]
def grad(g):
    dg = np.array([0,1])
    return dg

# Implementation of the crossing with the section with n_crossing = 1
    
lim = 2 * np.pi  # Ajuste a 2 pi

# Crear una malla de puntos usando meshgrid
x = np.linspace(-lim, lim, 100)  # 100 puntos en el eje x
y = np.linspace(-lim, lim, 100)  # 100 puntos en el eje y
X0, Y0 = np.meshgrid(x, y)  # Generar la malla

# Vector de las velocidades, paso como argumento todos los puntos de mi cuadrícula
u, v = oscilador_harmonico(0, [X0, Y0])  

# Lineas de corriente de nuestro sistema, que forma tiene
plt.figure(figsize=(8, 8)) 
plt.streamplot(X0, Y0, u, v, color='blue') 
plt.title('Líneas de corriente del oscilador armónico')  
plt.xlabel('x')  
plt.ylabel('y') 
plt.axhline(0, color='black', linewidth=0.5, ls='--')  
plt.axvline(0, color='black', linewidth=0.5, ls='--')  
plt.grid()  
plt.show()  
#%%
def Poincar_Map(initial_point, ncrossings, idir, function_to_evaluate):

    crossings = []
    results = []
    absolute_t = 0
  
    
 
    # The process repeats as much as many ncrossing we have
    for j in range(ncrossings):
        # We first integrate until we cross y=0 (backwards or forwards)
        t = 1
        

        t_periodic = [0, idir*2*np.pi]
        # Remark that for t = 0 we are on y = 0. This is not a good initial point, we need to perturbate the system
        perturbation_time = 1E-14
        # the perturbation time goes backward if idir == -1 and forward if idi == +1
        sol = solve_ivp(function_to_evaluate, t_periodic, initial_point, t_eval=np.array([idir * perturbation_time]), rtol=3e-14, atol=1e-14)
        # solve_ivp returns a OdeResult with the times t, the y (positions for each time)
        # for the positions the first array [0] is the x-coordinate, [1] y-coordinate. Then we want the first point as the initial point
        
        # Definition of the new initial time, the same as sol.t (this las varaible returnes as an np.array)
        absolute_t += perturbation_time
        
        # Coordinates of the new initial point
        x = sol.y[0][0]
        y = sol.y[1][0]
       
        
        # Now we start from an initial point and we know our systems is 2pi periodic
        # Points between the range (0, 2pi) we are gonna calculare the function
        iterations = 500
        teval = np.linspace(0, idir * 2 * np.pi, iterations)
        # Resolvemos para 500 puntos del periodo 2pi, hacia adelante o hacia atrás
        # rtol and atol, value parameter of the tolerance
        sol = solve_ivp(function_to_evaluate, t_periodic, y0=[x, y], t_eval=teval, rtol=3e-14, atol=1e-14)
        # Now we want to evaluate when the solution changes sign, that means we have crossed the y = 0
        t_ini = 0
        crossed = False
        for i in range(iterations): 
            point = sol.y[1][i]
            if idir == 1:     
                
                if (y*point)<0:
                    crossed = True
                    # we have a change of sign, y = 0 crossed
                    x_before, y_before = sol.y[0][i-1], sol.y[1][i-1]
                    x_after, y_after = sol.y[0][i], sol.y[1][i]
                    t_before, t_after = sol.t[i-1], sol.t[i]
                    crossx, crossy, t_ini = x_before, y_before, t_before
                    break
            if idir == -1:
                 
                if (y*point)<0:
                    crossed = True
                    # we have a change of sign, y = 0 crossed
                    x_before, y_before = sol.y[0][i-1], sol.y[1][i-1]
                    x_after, y_after = sol.y[0][i], sol.y[1][i]
                    t_before, t_after = sol.t[i-1], sol.t[i]
                    crossx, crossy, t_ini = x_after, y_after, t_after
                    break
        
        crossings.append([crossx, crossy])
        # Now we apply Newton method from t_before, not the same as initial time, to compure tm+1 = tm -G(tm)/G'(tm)
        tol = 10e-14
        belongs_to = False
        while not belongs_to:
            # we are doing this for each number of crossing times
            
            sol = solve_ivp(function_to_evaluate, [0,2*np.pi], y0=crossings[j], t_eval = np.array([t]),rtol=3e-14, atol=1e-14)
            t1 = t - (g(sol.y)/ (np.dot(grad(sol.y), function_to_evaluate(0, sol.y))))[0]
            t1 = t1 % (2 * np.pi)
            if np.abs(t - t1) < tol:
                belongs_to = True
            
            t = t1
            
        final = solve_ivp(function_to_evaluate, [0,2*np.pi], y0 = crossings[j], t_eval=np.array([t]),rtol=3e-14, atol=1e-14)
        # final point
        x = final.y[0][0]
        y = final.y[1][0]
        initial_point = x,y
        absolute_t += t + t_ini
        
        
        results.append([x,y, absolute_t])
    # Nos devuelve una lista con los puntos de corte y el tiempo en el que lo cruza
    return results
initial = [1,0]
# Poincar_Map(initial_point, ncrossings, idir, function_to_evaluate):
print("Forward: ")
result_forward = (Poincar_Map(initial, 2, 1, oscilador_harmonico))
print(result_forward)
time_second_crossing = result_forward[1][2]
print(f"time forward = {time_second_crossing}")
# Calculamos la diferencia con 2pi
difference = np.abs(2 * np.pi - time_second_crossing)

# Imprimir con 16 dígitos de precisión
print(f"difference with 2pi = {difference:.16f}")

print("Backward: ")
result_backward = (Poincar_Map(initial, 2, -1, oscilador_harmonico))
print(result_backward)
time_second_crossing_back = result_backward[1][2]
print(f"time backward = {time_second_crossing_back}")
# Calculamos la diferencia con 2pi
difference_back = np.abs(-2 * np.pi - time_second_crossing_back)

# Imprimir con 16 dígitos de precisión
print(f"difference with 2pi = {difference_back:.16f}")
#%%
# Now we want to use a initial point (0,1)

initial = [0,1]
# Poincar_Map(initial_point, ncrossings, idir, function_to_evaluate):
print("Forward: ")
result_forward = (Poincar_Map(initial, 2, 1, oscilador_harmonico))
print(result_forward)
time_second_crossing = result_forward[1][2]
print(f"time forward = {time_second_crossing}")
# Calculamos la diferencia con 2pi
difference = np.abs(3/2*np.pi - time_second_crossing)

# Imprimir con 16 dígitos de precisión
# print(f"difference with 2pi = {difference:.16f}")

print("Backward: ")
result_backward = (Poincar_Map(initial, 2, -1, oscilador_harmonico))
print(result_backward)
time_second_crossing_back = result_backward[1][2]
print(f"time backward = {time_second_crossing_back}")
# Calculamos la diferencia con 2pi
difference_back = np.abs(-3/2*np.pi- time_second_crossing_back)

print(difference_back)
print(difference)

# # Imprimir con 16 dígitos de precisión
# print(f"difference with 2pi = {difference_back:.16f}")

#%%
# PART B
def linear_system(x, y, a, b, c, d):
    result =  np.array([a*x+b*y, c*x+d*y])
    return result

# We are asked to plot the phase space for different parameters a,b,c and d

#%%
# CASE 1:
import numpy as np
import matplotlib.pyplot as plt

# Parámetros del sistema
a, b, c, d = 2, -5, 1, -2

# Generación de los valores de los ejes con np.linspace
lim = 2 * np.pi
x = np.linspace(-lim, lim, 500)
y = np.linspace(-lim, lim, 500)

# Crear las mallas
X0, X1 = np.meshgrid(x, y)

# Vector field del sistema
u, v = linear_system(X0, X1, a, b, c, d)

# Graficar el campo vectorial
fig = plt.axes()
fig.streamplot(X0, X1, u, v, density=0.6, broken_streamlines=False)
fig.set_aspect('equal', 'box')
plt.ylabel("y")
plt.xlabel("x")
plt.show()

#%%
# CASE 2:
import numpy as np
import matplotlib.pyplot as plt

# Parámetros del sistema
a, b, c, d = 3,-2,4, -1

# Generación de los valores de los ejes con np.linspace
lim = 2 * np.pi
x = np.linspace(-lim, lim, 500)
y = np.linspace(-lim, lim, 500)

# Crear las mallas
X0, X1 = np.meshgrid(x, y)

# Vector field del sistema
u, v = linear_system(X0, X1, a, b, c, d)

# Graficar el campo vectorial
fig = plt.axes()
fig.streamplot(X0, X1, u, v, density=0.6, broken_streamlines=False)
fig.set_aspect('equal', 'box')
plt.ylabel("y")
plt.xlabel("x")
plt.show()
#%%
# CASE 3:
import numpy as np
import matplotlib.pyplot as plt

# Parámetros del sistema
a, b, c, d = -1, 0, 3, 2

# Generación de los valores de los ejes con np.linspace
lim = 2 * np.pi
x = np.linspace(-lim, lim, 500)
y = np.linspace(-lim, lim, 500)

# Crear las mallas
X0, X1 = np.meshgrid(x, y)

# Vector field del sistema
u, v = linear_system(X0, X1, a, b, c, d)

# Graficar el campo vectorial
fig = plt.axes()
fig.streamplot(X0, X1, u, v, density=0.9, broken_streamlines=False)
fig.set_aspect('equal', 'box')
plt.ylabel("y")
plt.xlabel("x")
plt.show()



#%%
# The next part is to compute the eigenvalue and eigenvectors of the Jacobian matrix at (0,0)
# We are also asked to compute the stable and unstable manifold of the center


matrix = np.array([[a, b], [c, d]])
print(np.linalg.eig(matrix))
a, b, c, d = -1, 0, 3, 2
# Our system
def linear_system2(t, x):
    a, b, c, d = -1, 0, 3, 2
    return linear_system(x[0], x[1], a, b, c, d)

# MANIFOLDS
X = [1, 1]
X[1], X[0] = np.mgrid[-6:6:200j, -6:6:200j]
# flow of our phase space
u, v = linear_system(X[0], X[1], a, b, c, d)
fig = plt.axes()
fig.streamplot(X[0], X[1], u, v, density=0.4, broken_streamlines=False, color = "royalblue")
fig.set_aspect('equal', 'box')

fig.scatter(0, 0, label = "equilibrium point", color = "black")
s= 10e-7
# Compute and plot Wu+
x0 = np.array([0, 0]) + s*np.linalg.eig(matrix)[1][:, 0]
sol = solve_ivp(linear_system2, [0, 100], x0, t_eval=np.linspace(0, 10, 200), rtol=3e-14, atol=1e-14)
fig.plot(sol.y[0], sol.y[1], color='red', label = "unstable")
# Compute and plot Wu
x0 = np.array([0, 0]) - s*np.linalg.eig(matrix)[1][:, 0]
sol = solve_ivp(linear_system2, [0, 100], x0, t_eval=np.linspace(0, 10, 200), rtol=3e-14, atol=1e-14)
fig.plot(sol.y[0], sol.y[1], color='red')
# Compute and plot Ws+
x0 = np.array([0, 0]) + s*np.linalg.eig(matrix)[1][:, 1]
sol = solve_ivp(linear_system2, t_span=(0, -20), y0=x0, t_eval=np.linspace(0, -20, 200), rtol=3e-14, atol=1e-14)
fig.plot(sol.y[0], sol.y[1], color='yellow', label = "stable")
# Compute and plot Ws
x0 = np.array([0, 0]) - s*np.linalg.eig(matrix)[1][:, 1]
sol = solve_ivp(linear_system2, t_span=(0, -20), y0=x0, t_eval=np.linspace(0, -20, 200), rtol=3e-14, atol=1e-14)
fig.plot(sol.y[0], sol.y[1], color='yellow')

plt.xlim(-6, 6)
plt.ylim(-6, 6)
plt.xlabel("x")
plt.ylabel("y")
plt.legend(loc = "upper left")

plt.show()



#%%
# PART C:PENDULUM
# We are asked to plot the phase portrait of a pendulum, the ORBITS and DIRECTIONS
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
# arbritary parameters, the simplest ones
m = 1
g = 1
l = 1
def pendulum(x):
    # x[0] es el ángulo theta (posición angular)
    theta = x[0]
    
    # x[1] es la velocidad angular omega
    omega = x[1]
    
    # dtheta_dt es la derivada de theta con respecto al tiempo, es decir, la velocidad angular (omega)
    dtheta_dt = omega
    
    # domega_dt es la derivada de omega con respecto al tiempo, es decir, la aceleración angular
    # Se calcula usando la ecuación del movimiento del péndulo: -g/l * sin(theta)
    domega_dt = -(g / l) * np.sin(theta)
    
    # Devolver las derivadas como un array, d/dt theta y d/dt omega
    return dtheta_dt, domega_dt



# Parámetros
lim = 2 * np.pi
initial_point = [0, 0]
X = initial_point

# all of initial points we are gonna use
X[1], X[0] = np.mgrid[-lim:lim:100j, -2*lim:2*lim:100j]

# derivatives for omega and theta
u, v = pendulum(X)


fig, ax = plt.subplots(figsize=(8, 6))

# we want the phase portrait and their diretions
ax.streamplot(X[0], X[1], u, v, density=0.6, broken_streamlines=False, color="green")


ax.set_aspect('equal', 'box')

x_vals = [0, np.pi, -np.pi, 2*np.pi, -2*np.pi]
y_vals = [0, 0, 0, 0, 0]


ax.scatter(x_vals, y_vals, color='orange', label=r'$k\pi$')


ax.legend(loc='upper left')

ax.set_xlim([-3/2*lim, 3/2*lim])
ax.set_ylim([-lim, lim])
def omega_integral(theta, H):
    omega_p, omega_n =  np.sqrt(2 * (H - m*g*l*(1 - np.cos(theta))) / (m*l**2)), -np.sqrt(2 * (H - m*g*l*(1 - np.cos(theta))) / (m*l**2))
    return omega_p, omega_n

# H = T+V
x = np.linspace(-4*np.pi, 4*np.pi, 5000)
# Energy to achieve theta = -pi, calculating with our parameter H = 2 = 2mgl
H = 2
y_p, y_n = omega_integral(x, H)
plt.plot(x,y_p, label = "positive", color = "red")
plt.plot(x,y_n, label = "negative", color = "royalblue")


# Mostrar el gráfico
plt.grid(True)






#%%
# PART D
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Tamaño del plot
size = 5  

# Crear una malla de puntos en un rango determinado
x_values, y_values = np.meshgrid(np.linspace(0, size, 100), np.linspace(0, size, 100))

# Definir el sistema de ecuaciones Lotka-Volterra
def lotka_volterra(t, z):
    x, y = z
    dxdt = x * (3 - x - 2 * y)
    dydt = y * (2 - x - y)
    return [dxdt, dydt]

# Evaluar el sistema para cada punto de la malla
def calcular_vector_field(x_values, y_values):
    u = np.zeros_like(x_values)
    v = np.zeros_like(y_values)
    
    for i in range(x_values.shape[0]):
        for j in range(x_values.shape[1]):
            u[i, j], v[i, j] = lotka_volterra(0, [x_values[i, j], y_values[i, j]])
    
    return u, v

# Calcular el campo vectorial en la malla
u, v = calcular_vector_field(x_values, y_values)

# Crear la figura y trazar el flujo usando streamplot
fig, ax = plt.subplots()
ax.streamplot(x_values, y_values, u, v, density=0.4, broken_streamlines=False)
ax.set_aspect('equal')

# Equilibrium points at the plot
equilibrium_points = np.array([[0, 0], [0, 2], [3, 0], [1, 1]])
scatter = ax.scatter(equilibrium_points[:, 0], equilibrium_points[:, 1], color = "red", label='Equilibrium Points')

# Agregar etiquetas para los ejes
ax.set_xlabel('x')
ax.set_ylabel('y')

plt.title('Phase Portrait', fontsize=12)

# Configurar límites de los ejes para que vayan de 0 a 5
plt.xlim(0, 5)
plt.ylim(0, 5)

# Añadir una leyenda que explique que los puntos rojos son los puntos de equilibrio


# Definir el Jacobiano
def jacobian_matrix(x):
    return np.array([[3 - 2*x[0] - 2*x[1], -2*x[0]], [-x[1], 2 - x[0] - 2*x[1]]])

# MANIFOLDS

# Wu (órbita inestable) - positivo
x0 = np.array([1, 1]) + 1E-7 * np.linalg.eig(jacobian_matrix(np.array([1, 1])))[1][:, 0]
sol = solve_ivp(lotka_volterra, [0, 100], x0, t_eval=np.linspace(0, 100, 200))
ax.plot(sol.y[0], sol.y[1], c='orange', label="Unstable Manifold")

# Wu (órbita inestable) - negativo
x0 = np.array([1, 1]) - 1E-7 * np.linalg.eig(jacobian_matrix(np.array([1, 1])))[1][:, 0]
sol = solve_ivp(lotka_volterra, [0, 100], x0, t_eval=np.linspace(0, 100, 200))
ax.plot(sol.y[0], sol.y[1], c='orange')

# Ws (órbita estable) - positivo
x0 = np.array([1, 1]) + 1E-7 * np.linalg.eig(jacobian_matrix(np.array([1, 1])))[1][:, 1]
sol = solve_ivp(lotka_volterra, t_span=(0, -20), y0=x0, t_eval=np.linspace(0, -10, 200))
ax.plot(sol.y[0], sol.y[1], c='green', label="Stable Manifold")

# Ws (órbita estable) - negativo
x0 = np.array([1, 1]) - 1E-7 * np.linalg.eig(jacobian_matrix(np.array([1, 1])))[1][:, 1]
sol = solve_ivp(lotka_volterra, t_span=(0, -20), y0=x0, t_eval=np.linspace(0, -10, 200))
ax.plot(sol.y[0], sol.y[1], c='green')
plt.legend(loc='upper right')
# Mostrar la gráfica
plt.show()
#%%






