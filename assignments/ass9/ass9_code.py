# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 10:31:22 2024

@author: polop
"""
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
mu_range = [0.011720,0.011730,0.011740]


# mu_range = [0.304]

# Inicializar listas para almacenar resultados
x_list = []
y_list = []
s = 1e-6 # for smaller values i get an error
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

# # Coordenadas del segundo primario P2 para mu = 0.0202
# mu_target = 0.0202
# P2_x = 1 - mu_target
# P2_y = 0  # El segundo primario está en y = 0

# # Dibujar el punto P2 en el gráfico
# plt.scatter(-P2_x, -P2_y, color='red', marker='o', s=100, label=r'$P_2$ for $\mu=0.0202$')
plt.legend()
plt.show()



#%%  

import numpy as np
import matplotlib.pyplot as plt
from rich.progress import Progress
from scipy.integrate import solve_ivp
import os      

# Definir intervalos con nombres modificados
range_A = np.linspace(0.001, 0.015, round((0.015 - 0.001) / 0.00001) + 1)
range_B = np.linspace(0.015, 0.05, round((0.05 - 0.015) / 0.0001) + 1)
range_C = np.linspace(0.05, 0.49, round((0.49 - 0.05) / 0.001) + 1)

# Función para calcular los cruces
def compute_intersections(mu_val):
    for idx in range(len(get_eigenvalues(mu_val))):
        if np.isreal(get_eigenvalues(mu_val)[idx]) and get_eigenvalues(mu_val)[idx] > 0:
            eig_vector = get_eigenvectors(mu_val)[:, idx]

    initial_state = np.array([L3(mu_val), 0, 0, 0]) + 1E-6 * eig_vector
    if initial_state[1] > 0:
        initial_state = np.array([L3(mu_val), 0, 0, 0]) - 1E-6 * eig_vector

    results = poincare_map(initial_state, 1, mu_val)[0]
    return np.array([results[0][0], results[1][0], results[2][0], results[3][0]])

# Inicializar listas para almacenar resultados
res_A_x, res_A_y, res_A_xp, res_A_yp = [], [], [], []
res_B_x, res_B_y, res_B_xp, res_B_yp = [], [], [], []
res_C_x, res_C_y, res_C_xp, res_C_yp = [], [], [], []

# Procesar intervalos con `rich.progress`
with Progress() as progress:
    task_A = progress.add_task("[blue]Processing Range A...", total=len(range_A))
    for mu_val in range_A:
        results = compute_intersections(mu_val)
        res_A_x.append(results[0])
        res_A_y.append(results[1])
        res_A_xp.append(results[2])
        res_A_yp.append(results[3])
        progress.update(task_A, advance=1)

    task_B = progress.add_task("[green]Processing Range B...", total=len(range_B))
    for mu_val in range_B:
        results = compute_intersections(mu_val)
        res_B_x.append(results[0])
        res_B_y.append(results[1])
        res_B_xp.append(results[2])
        res_B_yp.append(results[3])
        progress.update(task_B, advance=1)

    task_C = progress.add_task("[magenta]Processing Range C...", total=len(range_C))
    for mu_val in range_C:
        results = compute_intersections(mu_val)
        res_C_x.append(results[0])
        res_C_y.append(results[1])
        res_C_xp.append(results[2])
        res_C_yp.append(results[3])
        progress.update(task_C, advance=1)

# Función para guardar los resultados en un archivo .txt
def save_results_to_file(filename, data_ranges, results, labels):
    with open(filename, 'w') as output:
        for range_data, result_data, label in zip(data_ranges, results, labels):
            output.write(f"# Results for {label}\n")
            for val_mu, val_result in zip(range_data, result_data):
                output.write(f"{val_mu:.6f} {val_result:.6f}\n")
            output.write("\n")
    print(f"Data has been saved to {os.path.abspath(filename)}")

# Guardar los resultados para los tres intervalos
output_filename = "computed_intersections.txt"
data_intervals = [range_A, range_B, range_C]
results_list = [res_A_xp, res_B_xp, res_C_xp]
interval_labels = ["Interval A", "Interval B", "Interval C"]

save_results_to_file(output_filename, data_intervals, results_list, interval_labels)

# Graficar los resultados
plt.figure(figsize=(10, 6))
plt.scatter(range_B, res_B_xp, color="brown", label="Range B")
plt.scatter(range_C, res_C_xp, color="brown", label="Range C")
plt.scatter(range_A, res_A_xp, color="brown", label="Range A")
plt.xlabel(r"$\mu$", fontsize=14)
plt.ylabel(r"$x'$", fontsize=14)
plt.title("Combined Plot of $x'$ vs $\mu$")
plt.legend()
plt.grid(True)
plt.show()
