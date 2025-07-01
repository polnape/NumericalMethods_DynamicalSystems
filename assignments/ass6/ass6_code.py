# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 15:38:38 2024

@author: polop
"""

import numpy as np
from scipy.integrate import ode

# Constants
mu = 0.012
x0 = [1.119748482176185, 0, 0, -0.2260036463376957]  # Initial condition
dir = 1  # 1 for forwards, -1 for backwards in time
n_crossings = 2  # Number of intersection points (crossings)
h = 0.01  # Integration time
tol = 1e-13  # Tolerance
period = 6.192169331318165  # Example period (from previous assignments)

# Function f for the differential equations (RTBP)
def f(t, x, mu, dir):
    r1 = np.sqrt((x[0] - mu)**2 + x[1]**2)
    r2 = np.sqrt((x[0] - mu + 1)**2 + x[1]**2)
    
    df = np.zeros(4)
    df[0] = x[2]
    df[1] = x[3]
    df[2] = 2*x[3] + x[0] - ((1-mu)*(x[0]-mu)/(r1**3)) - mu*(x[0]-mu+1)/(r2**3)
    df[3] = -2*x[2] + x[1]*(1 - (1-mu)/(r1**3) - mu/(r2**3))
    
    if dir == -1:
        df = -df
    return df

# Function g(x) which defines Sigma (y=0 plane, so g(x) = y)
def g(x):
    return x[1]

# Newton's method to find a better approximation of the intersection
def DeltaT(x0, mu):
    return -g(x0) / np.dot([0, 1, 0, 0], f(0, x0, mu, 1))

# First integral (Jacobi constant) calculation
def jacobi_constant(x, mu):
    r1 = np.sqrt((x[0] - mu)**2 + x[1]**2)
    r2 = np.sqrt((x[0] - mu + 1)**2 + x[1]**2)
    return x[0]**2 + x[1]**2 + 2*(1-mu)/r1 + 2*mu/r2 - (x[2]**2 + x[3]**2)

# Set up Poincaré section computation
poinc = np.zeros((n_crossings, 4))  # Store intersection points
total_time = 0
xk = np.copy(x0)  # Initialize xk with the initial condition

jacobi_initial = jacobi_constant(x0, mu)  # Initial Jacobi constant

# ODE solver setup
for i in range(n_crossings):
    # Define x0 as the last intersection point found (except for the first iteration)
    if i != 0:
        x0 = xk

    # Initialize the ODE solver
    solver = ode(f).set_integrator('dopri5', atol=tol, rtol=tol)
    solver.set_initial_value(x0, 0).set_f_params(mu, dir)
    
    # Find a point close to Sigma by integrating h seconds each time
    bool_crossed = False
    while not bool_crossed:
        solver.integrate(solver.t + h)
        total_time += h
        xk1 = solver.y
        
        # Check if we've crossed Sigma
        if g(x0) * g(xk1) < 0:
            bool_crossed = True
        x0 = xk1

    # Refine the intersection using Newton's method
    while abs(g(x0)) > tol:
        deltaT = DeltaT(x0, mu)
        total_time += dir * deltaT
        solver.set_initial_value(x0, 0).set_f_params(mu, np.sign(deltaT))
        solver.integrate(solver.t + abs(deltaT))
        x0 = solver.y

    # Store the intersection point
    xk = np.copy(x0)
    poinc[i, :] = xk

# Check final Jacobi constant
jacobi_final = jacobi_constant(xk, mu)

# Comparison between initial and final points (12 digit precision)
initial_vs_final = np.allclose(x0, xk, atol=1e-12)

# Check time of integration (difference must be less than 10^(-11))
time_difference = np.abs(total_time - period) < 1e-11
print(time_difference)

# Output results
xk, poinc, jacobi_initial, jacobi_final, initial_vs_final, total_time, time_difference
