import numpy as np
import matplotlib.pyplot as plt
from weno_lib import *
# ---------------------------------------------------------
# 1. Setup and Equation of State (Euler Flow)
# ---------------------------------------------------------
gamma = 1.4
CFL = 0.8
t_end = 0.2
nx = 200
x = np.linspace(-0.5, 0.5, nx)
dx = x[1] - x[0]

def get_pressure(U):
    rho = U[:, 0]
    m = U[:, 1]
    E = U[:, 2]
    v = m / rho
    return (gamma - 1.0) * (E - 0.5 * rho * v**2)

def get_flux(U):
    rho = U[:, 0]
    m = U[:, 1]
    E = U[:, 2]
    v = m / rho
    p = get_pressure(U)
    
    F = np.zeros_like(U)
    F[:, 0] = m
    F[:, 1] = rho * v**2 + p
    F[:, 2] = v * (E + p)
    return F

def get_max_speed(U):
    rho = U[:, 0]
    m = U[:, 1]
    p = get_pressure(U)
    v = m / rho
    c = np.sqrt(gamma * p / rho)
    return np.abs(v) + c

def get_jacobian(U):
    """
    Computes the characteristic matrix A = dF/dU evaluated at zone centers.
    """
    rho = U[:, 0]
    m = U[:, 1]
    E = U[:, 2]
    v = m / rho
    p = get_pressure(U)
    H = (E + p) / rho
    
    A = np.zeros((U.shape[0], 3, 3))
    A[:, 0, 1] = 1.0
    A[:, 1, 0] = 0.5 * (gamma - 3.0) * v**2
    A[:, 1, 1] = (3.0 - gamma) * v
    A[:, 1, 2] = gamma - 1.0
    A[:, 2, 0] = v * (0.5 * (gamma - 1.0) * v**2 - H)
    A[:, 2, 1] = H - (gamma - 1.0) * v**2
    A[:, 2, 2] = gamma * v
    return A


# ---------------------------------------------------------
# 3. AFD-WENO Update Implementation (Section 5)
# ---------------------------------------------------------
def get_afd_weno_rhs(U):
    # Step (iii): Get boundary interpolated states and cell centered derivatives
    U_L_eval, U_R_eval, dUdx_center = weno_ao_5_3_interpolation(U, dx, 'edge')
    
    # Form left/right Riemann states at faces (i+1/2)
    U_L_face = U_R_eval[:-1]  # from cell i, evaluated at right edge
    U_R_face = U_L_eval[1:]   # from cell i+1, evaluated at left edge
    
    # Compute characteristic matrix A at cell centers
    A_centers = get_jacobian(U)
    
    # Compute f_i = (A * dUdx) at cell centers
    f_centers = np.einsum('nij,nj->ni', A_centers, dUdx_center)
    
    # Step (vii): Section 4.1 Zone Boundary Interpolation
    d1 = np.zeros((nx, 3))
    d3 = np.zeros((nx, 3))
    for comp in range(3):
        d1c, d3c = weno_ao_43_boundary(f_centers[:, comp], dx, boundary_mode = 'edge')
        d1[:, comp] = d1c
        d3[:, comp] = d3c
    # Step (v): Apply LLF Riemann solver
    F_L = get_flux(U_L_face)
    F_R = get_flux(U_R_face)
    s_L = get_max_speed(U_L_face)
    s_R = get_max_speed(U_R_face)
    alpha = np.maximum(s_L, s_R)[:, np.newaxis]
    
    F_LLF = 0.5 * (F_L + F_R - alpha * (U_R_face - U_L_face))
    
    F_num = F_LLF - (dx**2 / 24.0) * d1[:-1] + (7 * dx**4 / 5760) * d3[:-1]
    
    # Return -1/dx * (F_{i+1/2} - F_{i-1/2})
    dUdt = np.zeros_like(U)
    dUdt[1:-1] = -(F_num[1:] - F_num[:-1]) / dx
    return dUdt

# ---------------------------------------------------------
# 4. Main Simulation Loop (Sod Shock Tube - Test 1)
# ---------------------------------------------------------
# Initial conditions
U = np.zeros((nx, 3))
left_mask = x < 0
right_mask = x >= 0

# Left state (rho, v, p) = (1.0, 0.0, 1.0)
# U = (rho, rho*v, E)
U[left_mask, 0] = 1.0
U[left_mask, 1] = 0.0
U[left_mask, 2] = 1.0 / (gamma - 1.0) # E = p/(g-1) + 0.5*rho*v^2

# Right state (rho, v, p) = (0.125, 0.0, 0.1)
U[right_mask, 0] = 0.125
U[right_mask, 1] = 0.0
U[right_mask, 2] = 0.1 / (gamma - 1.0)

t = 0.0
while t < t_end:
    max_s = np.max(get_max_speed(U))
    dt = CFL * dx / max_s
    if t + dt > t_end:
        dt = t_end - t
        
    # SSP-RK3 Time Stepping
    U_1 = U + dt * get_afd_weno_rhs(U)
    U_2 = 0.75 * U + 0.25 * U_1 + 0.25 * dt * get_afd_weno_rhs(U_1)
    U = (1.0/3.0) * U + (2.0/3.0) * U_2 + (2.0/3.0) * dt * get_afd_weno_rhs(U_2)
    
    t += dt

# Plot results
p = get_pressure(U)
v = U[:, 1] / U[:, 0]
rho = U[:, 0]

plt.figure(figsize=(12, 4))
plt.subplot(1, 3, 1)
plt.plot(x, rho, 'o-', markersize=3, color='k')
plt.title('Density')

plt.subplot(1, 3, 2)
plt.plot(x, v, 'o-', markersize=3, color='k')
plt.title('Velocity')

plt.subplot(1, 3, 3)
plt.plot(x, p, 'o-', markersize=3, color='k')
plt.title('Pressure')
plt.tight_layout()
plt.show()
