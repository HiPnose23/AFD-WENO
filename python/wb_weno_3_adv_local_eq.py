import numpy as np
import matplotlib.pyplot as plt
from weno_lib import *

# =========================================================
# 1. Setup and Parameters
# =========================================================
c_speed = 1.0    # Advection speed
lam = -0.2       # Reaction rate
CFL = 0.8
t_end = 2.0
nx = 100
x, dx = np.linspace(-1.0, 1.0, nx, endpoint=False, retstep=True)

# Base Equilibrium Profile: E(x) = exp((lam/c) * x)
E_center = np.exp((lam / c_speed) * x)
x_face = x + dx / 2.0
E_face = np.exp((lam / c_speed) * x_face)
E_prime_center = (lam / c_speed) * E_center

# =========================================================
# 2. Well-Balanced Flux Function 
# =========================================================
def compute_wb_flux(u):
    """Computes the 3rd-order AFD-WENO flux using equilibrium variables."""
    # Transform to equilibrium variable w
    w = u / E_center
    
    # WENO-AO(3) interpolation on w
    w_L_eval, w_R_eval, dwdx_center = weno_ao_3_interpolation(w, dx, 'wrap')
    w_L_face = w_R_eval
    w_R_face = np.roll(w_L_eval, -1)
    
    # Transform back to original variable u
    u_L_face = w_L_face * E_face
    u_R_face = w_R_face * E_face
    dudx_center = dwdx_center * E_center + w * E_prime_center
    
    # Upwind Riemann Flux
    F_star = c_speed * u_L_face if c_speed >= 0 else c_speed * u_R_face
        
    # f = A * dudx
    f_centers = c_speed * dudx_center
    
    # Boundary derivative of f
    d1dx, _ = weno_ao_43_boundary(f_centers, dx, boundary_mode='wrap')
    
    # 3rd-order AFD-WENO flux correction
    F_num = F_star - (dx**2 / 24.0) * d1dx
    return F_num

# =========================================================
# 3. Local Equilibrium RHS 
# =========================================================
def get_rhs_local_equilibrium(u):
    """
    Evaluates the source term cell-by-cell by constructing a local
    equilibrium state and passing it through the numerical flux operator.
    """
    # 1. Compute physical flux and gradient for the actual state
    F_num_actual = compute_wb_flux(u)
    dFdx_actual = (F_num_actual - np.roll(F_num_actual, 1)) / dx
    
    S_discrete = np.zeros(nx)
    
    # 2. Cell-by-cell Local Equilibrium Source Term
    for i in range(nx):
        # Step A: Find the local scalar for cell i
        w_local = u[i] / E_center[i]
        
        # Step B: Build the "imaginary perfect world" matching cell i
        u_eq_local = w_local * E_center
        
        # Step C: Pass this imaginary world through the exact flux solver
        F_num_imaginary = compute_wb_flux(u_eq_local)
        
        # Step D: Define the source term for cell i as the flux gradient 
        # of this imaginary world exactly across cell i.
        S_discrete[i] = (F_num_imaginary[i] - F_num_imaginary[i-1]) / dx
        
    # 3. Final RHS
    dudt = -dFdx_actual + S_discrete
    return dudt

# =========================================================
# 4. Verification Test
# =========================================================
print("\n--- EXACT STEADY-STATE PRESERVATION TEST ---")
print("Approach: Local Equilibrium Evaluation (General Xing-Shu)")

# Initialize with exact steady state
u = E_center.copy() 

# Evaluate RHS
rhs = get_rhs_local_equilibrium(u)

# Measure error in the interior
interior = slice(5, -5)
max_error = np.max(np.abs(rhs[interior]))

print(f"Max RHS Error in interior at t=0: {max_error:.3e}")

if max_error < 1e-12:
    print("SUCCESS: The scheme is well-balanced to machine precision!")
    print("The imaginary fluxes perfectly matched the actual fluxes.")
else:
    print("FAILED: Truncation error is still present.")

# ---------------------------------------------------------
# 3. Main Simulation Loop
# ---------------------------------------------------------
# Set initial condition exactly to the steady-state
u = E_center.copy() 
u_init = u.copy()

initial_dudt = get_rhs_local_equilibrium(u)
print(f"Max RHS error in interior at t=0: {np.max(np.abs(initial_dudt[1:])):.3e}")
# ---------------------------------

t = 0.0
while t < t_end:
    dt = CFL * dx / np.abs(c_speed)
    if t + dt > t_end:
        dt = t_end - t
        
    # SSP-RK3 Time Stepping
    u1 = u + dt * get_rhs_local_equilibrium(u)
    u2 = 0.75 * u + 0.25 * u1 + 0.25 * dt * get_rhs_local_equilibrium(u1)
    u = (1.0/3.0) * u + (2.0/3.0) * u2 + (2.0/3.0) * dt * get_rhs_local_equilibrium(u2)    
    t += dt


plt.figure(figsize=(8, 5))
plt.plot(x, u_init, '--', label='Exact Equilibrium (t=0)', color='gray', linewidth=2)
plt.plot(x, u, '-', label=f'WB-AFD-WENO (t={t_end})', color='blue')
plt.title(f'Well-Balanced Advection-Reaction ($\\lambda$={lam})')
plt.xlabel('x')
plt.ylabel('u')
plt.legend()
plt.grid(True)
plt.show()
