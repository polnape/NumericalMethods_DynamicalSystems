# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 14:40:27 2024

@author: polop
"""

# BISECTION METHOD TO FIND ROOTS

"""Find two points, say a and b such that a < b and f(a)* f(b) < 0
Find the midpoint of a and b, say “t”
t is the root of the given function if f(t) = 0; else follow the next step
Divide the interval [a, b] – If f(t)*f(a) <0, there exist a root between t and a
– else if f(t) *f (b) < 0, there exist a root between t and b
Repeat above three steps until f(t) = 0"""


import numpy as np
import matplotlib.pyplot as plt
from rich.progress import Progress
def bisection_method(f, x1, x2, tol=1e-8, max_iter=1000):

    # Caso especial: la raíz es exactamente uno de los límites
    if f(x1) == 0:
        return x1
    if f(x2) == 0:
        return x2

    # Comprobar condición de cambio de signo
    if f(x1) * f(x2) > 0:
        raise ValueError("El método de bisección requiere que f(x1) * f(x2) < 0.")

    iter_count = 0
    while abs(x2 - x1) > tol and iter_count < max_iter:
        c = (x1 + x2) / 2  # Punto medio del intervalo
        if f(c) == 0 or abs(x2 - x1) < tol:  # Raíz encontrada o intervalo suficientemente pequeño
            return c
        elif f(x1) * f(c) < 0:
            # raiz en la parte izquierda
            x2 = c  # La raíz está en [x1, c]
        else:
            # raiz en la parte derecha
            x1 = c  # La raíz está en [c, x2]
        # print(f(c))
        iter_count += 1

    return (x1 + x2) / 2  # Aproximación de la raíz después de alcanzar la tolerancia


# pip install rich necessary
from rich.progress import Progress
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import os




# Computation of Li points, using an iteration method for the quintic polynomial equation



def L1(mu):  # Returns the x for L1 given mu
    convergence = 0
    # initial value
    x = (mu / (3 * (1 - mu))) ** (1 / 3)
    tope_iter = 1000
    iteration_count = 0
    # itera method inside the own computation of the equilibrium point
    while convergence == 0 and iteration_count < tope_iter:
        xk = (mu * (1 - x) ** 2 / (3 - 2 * mu - x * (3 - mu - x))) ** (1 / 3)
        if np.abs(x - xk) < 1e-17:
            convergence = 1
        x = xk
        iteration_count += 1
    return mu - 1 + x

def L2(mu):  # Returns the x for L2 given mu
    convergence = 0    
    tope_iter = 1000
    iteration_count = 0
    # initial value epsilon_0
    x = (mu / (3 * (1 - mu))) ** (1 / 3)
    while convergence == 0 and iteration_count < tope_iter:
        xk = (mu * (1 + x) ** 2 / (3 - 2 * mu + x * (3 - mu + x))) ** (1 / 3)
        if np.abs(x - xk) < 1e-17:
            convergence = 1
        x = xk
        iteration_count += 1
    return mu - 1 - x

def L3(mu):  # Returns the x for L3 given mu
    convergence = 0
    x = 1 - 6 * mu / 12
    tope_iter = 1000
    iteration_count = 0
    while convergence == 0 and iteration_count < tope_iter:
        xk = ((1 - mu) * (1 + x) ** 2 / (1 + 2 * mu + x * (2 + mu + x))) ** (1 / 3)
        if np.abs(x - xk) < 1e-17:
            convergence = 1
        x = xk
        iteration_count += 1
    return mu + x


# Jacobi constant and derivatives:


# PARTIAL DERIVATIVES
def Omega_x(x, y, mu):
    r1 = np.sqrt((x - mu) ** 2 + y ** 2)
    r2 = np.sqrt((x - mu + 1) ** 2 + y ** 2)
    return x - (1 - mu) * (x - mu) / (r1 ** 3) - mu * (x - mu + 1) / (r2 ** 3)

def Omega_y(x, y, mu):
    r1 = np.sqrt((x - mu) ** 2 + y ** 2)
    r2 = np.sqrt((x - mu + 1) ** 2 + y ** 2)
    return y - (1 - mu) * y / (r1 ** 3) - mu * y / (r2 ** 3)

def Omega_xx(x, y, mu):
    r1 = np.sqrt((x - mu) ** 2 + y ** 2)
    r2 = np.sqrt((x - mu + 1) ** 2 + y ** 2)
    return 1 - (1 - mu) / r1 ** 3 + (3 * (1 - mu) * (x - mu) ** 2) / (r1 ** 5) - mu / r2 ** 3 + (3 * mu * (x - mu + 1) ** 2) / r2 ** 5

def Omega_xy(x, y, mu):
    r1 = np.sqrt((x - mu) ** 2 + y ** 2)
    r2 = np.sqrt((x - mu + 1) ** 2 + y ** 2)
    return y * ((3 * (1 - mu) * (x - mu)) / r1 ** 5 + (3 * mu * (x - mu + 1)) / r2 ** 5)

def Omega_yy(x, y, mu):
    r1 = np.sqrt((x - mu) ** 2 + y ** 2)
    r2 = np.sqrt((x - mu + 1) ** 2 + y ** 2)
    return 1 - (1 - mu) / r1 ** 3 - mu / r2 ** 3 + y ** 2 * ((3 * (1 - mu)) / r1 ** 5 + (3 * mu) / r2 ** 5)

# OMEGA TOTAL
def Omega(x, y, mu):
    r1 = np.sqrt((x - mu) ** 2 + y ** 2)
    r2 = np.sqrt((x - mu + 1) ** 2 + y ** 2)
    return 1 / 2 * (x ** 2 + y ** 2) + (1 - mu) / r1 + mu / r2 + 1 / 2 * mu * (1 - mu)

# Function to compute the dynamical_system


def dynamical_system(t, x, mu):
    # Cálculo de las distancias r1 y r2
    r1 = np.sqrt((x[0] - mu) ** 2 + x[1] ** 2)
    r2 = np.sqrt((x[0] - mu + 1) ** 2 + x[1] ** 2)
    
    # Cálculo de Omega_x y Omega_y
    Omega_x = x[0] - (1 - mu) * (x[0] - mu) / (r1 ** 3) - mu * (x[0] - mu + 1) / (r2 ** 3)
    Omega_y = x[1] - (1 - mu) * x[1] / (r1 ** 3) - mu * x[1] / (r2 ** 3)
    
    # Cálculo de las derivadas
    df1 = x[2]
    df2 = x[3]
    df3 = 2 * x[3] + Omega_x
    df4 = -2 * x[2] + Omega_y
    
    return df1, df2, df3, df4

# Jacobian of the RTBP
def Jacobian(x, y, mu):
    return np.array([
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [Omega_xx(x, y, mu), Omega_xy(x, y, mu), 0, 2],
        [Omega_xy(x, y, mu), Omega_yy(x, y, mu), -2, 0]
    ])

# Eigenvalues and eigenvectors of the L3
# remove [0] to get imaginary part 
def get_eigenvalues(mu):
    return np.linalg.eig(Jacobian(L3(mu), 0, mu))[0]

def get_eigenvectors(mu):
    return np.linalg.eig(Jacobian(L3(mu), 0, mu))[1]

# Functions for the Poincaré Map
def g_sigma(x):
    return x[1]

def grad_g_sigma(x):
    return np.array([0, 1, 0, 0])

import numpy as np
from scipy.integrate import solve_ivp

# Given an initial condition compute (x, 0, 0, yp)
def y_p(mu, x, C):
    # print(mu, x, C)
    y_prima = -1*np.sqrt(2*Omega(x, 0, mu)- C)
    
    x_p = 0
    y = 0 
    x = x
    state = [x, y, x_p, y_prima]
    return state

# Given (x, 0, 0, yp) computes f(x) = x'
def poincare_map(state, num_crossings, mu):
    
    # Determinar el límite de tiempo basado en el valor de mu
    def get_max_time(mu):
        ranges = {
            0.005: 500,
            0.01: 250,
            0.1: 100
        }
        for threshold, time in ranges.items():
            if mu < threshold:
                return time
        return 50

    max_time = get_max_time(mu)
    # max_time = 4

    # Número de puntos en el tiempo para evaluar
    num_eval_points = 500000

    # Comprobar que la dirección es válida
    direction_check = 1
    if direction_check not in [1, -1]:
        print("-1 BACKWARD, +1 FORWARD")
        return

    # Inicializar variables del sistema
    accumulated_time = 0
    pos_y = state[1]
    vel_x = state[2]
    vel_y = state[3]
    pos_x = state[0]

    # Iterar para buscar múltiples cruces
    for crossing in range(num_crossings):
        # Definir el intervalo de tiempo
        integration_span = [0, max_time]

        # Realizar un paso de integración corto si no es el primer cruce
        if crossing != 0:
            solution = solve_ivp(
                fun=lambda t, var: dynamical_system(t, var, mu),
                t_span=integration_span,
                y0=[pos_x, pos_y, vel_x, vel_y],
                t_eval=np.array([direction_check * 1E-10]),
                rtol=3e-14,
                atol=1e-14
            )
            accumulated_time += direction_check * 1E-10
            pos_x, pos_y, vel_x, vel_y = solution.y[:, 0]

        # Realizar la integración en un rango de tiempo
        evaluation_times = np.linspace(0, max_time, num_eval_points)
        solution = solve_ivp(
            fun=lambda t, var: dynamical_system(t, var, mu),
            t_span=integration_span,
            y0=[pos_x, pos_y, vel_x, vel_y],
            t_eval=evaluation_times,
            rtol=3e-14,
            atol=1e-14
        )

        initial_time = 0
        crossing_found = False

        # Detectar cruce con el eje y=0 hacia adelante o atrás
        if pos_y >= 0:
            for time_step in range(num_eval_points):
                if solution.y[1][time_step] < 0:
                    if direction_check == 1:
                        pos_x, pos_y, vel_x, vel_y = solution.y[:, time_step - 1]
                        initial_time = solution.t[time_step - 1]
                        crossing_found = True
                        break
                    else:
                        pos_x, pos_y, vel_x, vel_y = solution.y[:, time_step]
                        initial_time = solution.t[time_step]
                        crossing_found = True
                        break

        if pos_y < 0 and not crossing_found:
            for time_step in range(num_eval_points):
                if solution.y[1][time_step] > 0:
                    if direction_check == 1:
                        pos_x, pos_y, vel_x, vel_y = solution.y[:, time_step - 1]
                        initial_time = solution.t[time_step - 1]
                        break
                    else:
                        pos_x, pos_y, vel_x, vel_y = solution.y[:, time_step]
                        initial_time = solution.t[time_step]
                        break

        # Ajuste de Newton para refinar el cruce
        convergence_reached = False
        iterations = 0
        current_step_time = 0

        while not convergence_reached and iterations < 1000:
            solution = solve_ivp(
                fun=lambda t, var: dynamical_system(t, var, mu),
                t_span=integration_span,
                y0=[pos_x, pos_y, vel_x, vel_y],
                t_eval=np.array([current_step_time]),
                rtol=3e-14,
                atol=1e-14
            )

            t_adjusted = current_step_time - (g_sigma(solution.y) /
                                              (np.dot(grad_g_sigma(solution.y), dynamical_system(0, solution.y, mu))))[0]
            t_adjusted = np.real(t_adjusted) % (2 * np.pi)

            if np.abs(current_step_time - t_adjusted) < 1e-14:
                convergence_reached = True

            current_step_time = t_adjusted
            iterations += 1

            if iterations > 2000:
                print('Error: Too many iterations')

        solution = solve_ivp(
            fun=lambda t, var: dynamical_system(t, var, mu),
            t_span=integration_span,
            y0=[pos_x, pos_y, vel_x, vel_y],
            t_eval=np.array([current_step_time]),
            rtol=3e-14,
            atol=1e-14
        )

        pos_x, pos_y, vel_x, vel_y = solution.y[:, 0]
        accumulated_time += current_step_time + initial_time

    return [solution.y, accumulated_time]

# Definir el rango de mu
# mu_range = [0.005630, 0.005640, 0.005650]

# mu_range = [0.020000 
# ,0.020100 
# ,0.020200 ]
# mu_range = [0.009250, 
# 0.009260 ,
# 0.009270 ]


mu_range = [0.1]


# mu_range = [0.304]

# Inicializar listas para almacenar resultados
x_list = []
y_list = []
s = 1e-6 # for smaller values i get an error
i = 0
for mu in mu_range:
    # Encontrar el vector propio para el valor propio positivo
    for i in range(len(get_eigenvalues(mu))):
        if np.isreal(get_eigenvalues(mu)[i]) and get_eigenvalues(mu)[i] > 0:
            v0 = get_eigenvectors(mu)[:, i]
    
    # Calcular el punto inicial x0
    L3_value = L3(mu)
    x0 = np.array([L3_value, 0, 0, 0]) +  s* v0
    if x0[1] > 0:
        x0 = np.array([L3_value, 0, 0, 0]) - s * v0
    
    # Calcular el mapa de Poincaré y la solución del sistema
    a = poincare_map(x0, 1, mu)
    sol = solve_ivp(
        fun=lambda t, y: dynamical_system(t, y, mu),
        t_span=[0, a[1]],
        y0=x0,
        t_eval=np.linspace(0, a[1], 1000),
        rtol=3e-14,
        atol=1e-14
    )
    
    # Almacenar los resultados en las listas
    x_list.append(sol.y[0])
    y_list.append(sol.y[1])
   

# Graficar los resultados en un solo plot
plt.figure(figsize=(10, 6))
for j in range(len(mu_range)):
    plt.plot(x_list[j], y_list[j], label=f'$\mu$={mu_range[j]}')

plt.xlabel(r"$x$",fontsize = 14)
plt.ylabel(r"$y$",fontsize = 14)
# plt.xlim(-1.15, -0.95)
# plt.ylim(-0.04, 0.02)
plt.title("Poincaré Map Trajectories for different $\mu$ values")

plt.grid(True)
# Calcular el punto L2 para mu = 0.0202

# Coordenadas del segundo primario P2 para mu = 0.0202
mu_target = 0.0202
P2_x = 1 - mu_target
P2_y = 0  # El segundo primario está en y = 0

# Dibujar el punto P2 en el gráfico
plt.scatter(-P2_x, -P2_y, color='red', marker='o', s=100, label=r'$P_2$ for $\mu=0.0202$')
plt.legend()
plt.show() 

# Using the poincaré map we want to compute x' for a given x
def F(mu, x, C):
    # x as the initial condition
    #first find the yp to pass as the initial state 
    state = y_p(mu, x, C)
    solution = poincare_map(state, 2, mu)
    # sol.y y de aqui quiero la tercera columna x' y el primer valor incial de x'
    x_prima = solution[0][2][0]
    
    return x_prima

# PQ ES X' = 0 PARA N = 1???
C1 = 2.1
# C2 = 3.189
C2 = 2.99
mu = 0.1
L3_point = L3(mu)
print(L3_point)


#%%
# initial bisection points
x2 = 1.393505608247538
x1 = x2
fx1 = F(mu, x1, C2)
fx2 = F(mu, x2, C2)
print(x1, fx1)

# he inicializado los dos puntos en el mismo sitio, voy a desplazar uno a la derecha para ir iterando con mi método
# Limitamos el número de iteraciones para evitar bucles infinitos
max_iterations = 10000
for _ in range(max_iterations):
    # busco un cambio de signo en F(x)
    if (fx2 * fx1 <= 0):  # Condición de salida
        break
    x1 = x2
    x2 = x2+ 1e-3
    fx1 = fx2
    fx2 = F(0.1, x2, C2)
    print(f"x1 = {x1}, x2 = {x2}, fx1 = {fx1}, fx2 = {fx2}")
else:
    print("Error: No se encontró un cambio de signo en el intervalo dado.")

#%%
root = bisection_method(lambda x: F(0.1, x, C2), x1, x2)
print(root, F(0.1, root, C2))

# 1.0610139281023079 -4.119488789151957e-07

#%%
# # lo pongo para no volvero a correr
# root = 1.0610139281023079
root=  1.393505608247538
from rich.progress import Progress
# ITERATION FOR MULTPLE C, DELTA_C = 1E-4
num_points = int((C2 - C1) / 1e-3) + 1

# Crear un array equiespaciado con linspace
list_C = np.linspace(C2, C1, num_points)
mu = 0.1
# Cuantas orbitas habrá
root_orbits = np.zeros(len(list_C))
root_orbits[0] = root
# Proceso principal con medidor de progreso
with Progress() as progress:
    task = progress.add_task("[cyan]Processing...", total=len(list_C)-1)

    for i in range(1, len(list_C)):
        C_actual = list_C[i]
        x1 = root_orbits[i-1]
        x2 = x1
        fx1 = F(mu, x1, C_actual)
        fx2 = fx1
        
        while fx2*fx1>0:
            x2 = x1
            fx2 = fx1
            # Para C menores me desplazo a x más grandes
            x1 += 1e-3
            fx1 = F(mu, x1, C_actual)
            
        print(f"INTERVALO CONSEGUIDO: C_actual: {C_actual}, fx2:{ fx2}, fx1: {fx1}")
        nueva_raiz = bisection_method(lambda x: F(mu, x, C_actual), x1, x2)
        root_orbits[i] = nueva_raiz
            

        progress.update(task, advance=1)
    

#%%
with open("./results2.txt", "w") as file:
    file.write("roots_orbits\tc_list\n")  # Encabezado
    for root, c in zip(root_orbits, list_C):
        if root != 0:  # Solo guardar si el valor de root_orbits es distinto de 0
            file.write(f"{root}\t{c}\n")
     #%%
# Guardar list_C y root_orbits en un archivo .txt
with open("results2.txt", "w") as file:
    file.write("C values and their corresponding root orbits:\n")
    file.write("C_value\tRoot_orbit\n")
    for C, root in zip(list_C, root_orbits):
        file.write(f"{C:.6f}\t{root:.6f}\n")

print("Results saved to 'results2.txt'")     
import os

print(f"File saved at: {os.path.abspath('results2.txt')}")
  #%%
  
import os
import numpy as np

# Cambiar al directorio donde se encuentra el script
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Inicializar listas para almacenar los datos combinados
C_total = []
roots_total = []

# Función para leer un archivo y procesar los datos
def read_file(file_path, reverse_columns=False, skip_header=False):
    C_values = []
    roots_values = []
    with open(file_path, "r") as file:
        if skip_header:
            next(file)  # Saltar la primera línea si es un encabezado
        for line in file:
            try:
                # Leer y convertir valores dependiendo del orden
                val1, val2 = map(float, line.split())
                if reverse_columns:
                    C_values.append(val1)
                    roots_values.append(val2)
                else:
                    roots_values.append(val1)
                    C_values.append(val2)
            except ValueError:
                # Saltar líneas mal formateadas
                print(f"Línea inválida en {file_path}: {line.strip()}")
    return C_values, roots_values

# Leer el primer archivo (saltando el encabezado)
C1, roots1 = read_file("results.txt", reverse_columns=False, skip_header=False)

# Leer el segundo archivo (columnas invertidas, sin encabezado)
C2, roots2 = read_file("results2.txt", reverse_columns=True)

# Combinar datos
C_total = np.array(C1 + C2)
roots_total = np.array(roots1 + roots2)

# # Imprimir resultados para verificar
# print("C_total:", C_total)
# print("roots_total:", roots_total)
#%%
import matplotlib.pyplot as plt

# Crear el gráfico C(roots)
plt.figure(figsize=(8, 6))
plt.plot(roots_total, C_total, linestyle='-', color='royalblue', label='C(roots)')

# Personalizar el gráfico
# plt.title("Plot of C(initial points x)")
plt.xlabel("x", fontsize = 16)
plt.ylabel("C", fontsize = 16)
plt.grid(True)
plt.legend()
plt.tight_layout()

# Mostrar el gráfico
plt.show()
#%%
import numpy as np

# Inicializar arrays
root_orbits = np.zeros(1000)  # Ajusta el tamaño según sea necesario
list_C = np.zeros(1000)

# Leer el archivo para cargar valores procesados si existe
try:
    processed_roots = []
    processed_c_list = []
    
    with open("results.txt", "r") as file:
        next(file)  # Saltar el encabezado
        for line in file:
            # Manejo de errores por seguridad
            try:
                root, c = map(float, line.split())
                processed_roots.append(root)
                processed_c_list.append(c)
            except ValueError:
                print(f"Error procesando la línea: {line.strip()}")

    # Convertir a numpy arrays y actualizar root_orbits y list_C
    num_processed = len(processed_roots)
    root_orbits[:num_processed] = np.array(processed_roots)
    list_C[:num_processed] = np.array(processed_c_list)

except FileNotFoundError:
    num_processed = 0  # Si el archivo no existe, empieza desde el principio

#%%
# Proceso principal con medidor de progreso
with Progress() as progress:
    task = progress.add_task("[cyan]Processing...", total=len(list_C) - num_processed)

    for i in range(num_processed, len(list_C)):
        C_actual = list_C[i]
        x1 = root_orbits[i - 1]
        x2 = x1
        fx1 = F(mu, x1, C_actual)
        fx2 = fx1
        
        # Encontrar intervalo para la raíz
        while fx2 * fx1 > 0:
            x2 = x1
            fx2 = fx1
            x1 += 1e-3
            fx1 = F(mu, x1, C_actual)
        
        print(f"INTERVALO CONSEGUIDO: C_actual: {C_actual}, fx2: {fx2}, fx1: {fx1}")
        nueva_raiz = bisection_method(lambda x: F(mu, x, C_actual), x1, x2)
        root_orbits[i] = nueva_raiz

        # Guardar en el archivo en cada iteración
        with open("./results.txt", "a") as file:
            file.write(f"{nueva_raiz}\t{C_actual}\n")
        
        progress.update(task, advance=1)


#%%
C =[3.15, 2.5, 2.1]

# Encontramos los índices de los valores de C en C_total
indices = np.where(np.isin(C_total, C))[0]
x_list = []
y_list = []
# Creamos la lista x_iniciales con los valores correspondientes en root_orbits
x_iniciales = roots_total[indices]
mu = 0.1
i = 0
for c_values in C:
    x0 = x_iniciales[i]
    # print(x0)
    state = y_p(mu, x0, c_values)
    val = poincare_map(state, 2, mu)[1]*2
    print(f"val: {val}")
    sol = solve_ivp(
        fun=lambda t, y: dynamical_system(t, y, mu),
        t_span=[0, val],
        y0=state,
        t_eval=np.linspace(0, val, 1000),
        rtol=3e-14,
        atol=1e-14
        )
    x_list.append(sol.y[0])
    y_list.append(sol.y[1])
    i = i+1
#%%   
# Graficar los resultados en un solo plot
plt.figure(figsize=(10, 6))
for j in range(len(C)):
    plt.plot(x_list[j], y_list[j], label = f"C = {C[j]}")

plt.xlabel(r"$x$",fontsize = 14)
plt.ylabel(r"$y$",fontsize = 14)
# plt.xlim(-1.15, -0.95)
# plt.ylim(-0.04, 0.02)
plt.title("Poincaré Map Trajectories for different $C$ values")

plt.grid(True)

plt.scatter(-(1-0.1), 0, color = "gray", label = "P2", s = 80)
plt.scatter(0.1, 0, color = "brown", label ="P1", s =80)
# Coordenadas del segundo primario P2 para mu = 0.0202
mu_target = 0.1
# P2_x = 1 - mu_target
# P2_y = 0  # El segundo primario está en y = 0

# Dibujar el punto P2 en el gráfico
plt.scatter(L3_value, 0, color='red', marker='o', s=100, label=r'$L_3$')
plt.legend()
plt.show() 
#%%
import os
import numpy as np
import matplotlib.pyplot as plt

# Cambiar automáticamente al directorio del script
script_dir = os.path.dirname(os.path.abspath(__file__))  # Ruta del script
os.chdir(script_dir)  # Cambiar directorio de trabajo al del script
print(f"Directorio cambiado a: {os.getcwd()}")

# Nombre del archivo .txt
txt_file = "variables_datos.txt"

try:
    # Leer datos desde el archivo .txt, separador = tabulador (\t)
    data = np.loadtxt(txt_file, delimiter='\t')
    
    # Dividir los datos en las columnas correspondientes
    C_t = data[:, 0]  # Primera columna
    x_t = data[:, 1]  # Segunda columna
    T = data[:, 2]    # Tercera columna
    
    print(f"Archivo '{txt_file}' leído correctamente y procesado.")
except Exception as e:
    print(f"Error al leer el archivo: {e}")
    C_t, x_t, T = None, None, None

# Generar el gráfico de T en función de x_t
if C_t is not None and x_t is not None and T is not None:
    plt.figure(figsize=(8, 6))
    plt.plot(x_t, T,  linestyle='--', label='T(x)', color = "royalblue")
    plt.xlabel('x', fontsize = 16)
    plt.ylabel('T', fontsize = 16)
    # plt.title('Gráfico de T en función de x_t')
    plt.grid(True)
    # plt.legend()

plt.axhline(y = 5.83, color = "brown", label = "T = 5.83", linestyle = "--") 
plt.legend()   
plt.show()




