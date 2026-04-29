import numpy as np
import matplotlib.pyplot as plt
from weno_lib import *

# ---------------------------------------------------------
# 1. Setup and Parameters
# ---------------------------------------------------------
c_speed = 1.0    # Advection speed
lam = -0.2       # Reaction rate
CFL = 0.8
t_end = 2.0      
nx = 200
x = np.linspace(-1.0, 1.0, nx, endpoint=False)
dx = x[1] - x[0]

# ---------------------------------------------------------
# 2. Well-Balanced Components (Xing-Shu Framework)
# ---------------------------------------------------------
# Exact equilibrium state: E(x) = exp((lam/c) * x)
# We evaluate this at cell centers and cell faces.
E_center = np.exp((lam / c_speed) * x)
x_face = x + dx / 2.0
E_face = np.exp((lam / c_speed) * x_face)

# Analytical derivative of the equilibrium state at cell centers
E_prime_center = (lam / c_speed) * E_center

def compute_wb_flux(u):
    """Computes the strictly 3rd-order AFD-WENO flux using equilibrium variables."""
    # 1. Transform to equilibrium variable w (which is a constant at steady state)
    w = u / E_center
    
    # 2. Interpolate w using WENO-AO(3)
    w_L_eval, w_R_eval, dwdx_center = weno_ao_3_interpolation(w, dx, 'wrap')
    
    # Boundary states at faces i+1/2
    w_L_face = w_R_eval
    w_R_face = np.roll(w_L_eval, -1)
    
    # 3. Transform back to original variable u
    u_L_face = w_L_face * E_face
    u_R_face = w_R_face * E_face
    
    # Product rule for the derivative: d(w*E)/dx = w'*E + w*E'
    dudx_center = dwdx_center * E_center + w * E_prime_center
    
    # 4. Upwind Riemann Flux at face i+1/2
    if c_speed >= 0:
        F_star = c_speed * u_L_face
    else:
        F_star = c_speed * u_R_face
        
    # 5. f = A * dudx (where A = c_speed)
    f_centers = c_speed * dudx_center
    
    # 6. Boundary derivative of f (strictly 3rd order, so we only need d1dx)
    d1dx, _ = weno_ao_43_boundary(f_centers, dx, boundary_mode='wrap')
    
    # 7. 3rd-order AFD-WENO high-order flux correction
    F_num = F_star - (dx**2 / 24.0) * d1dx
    return F_num

# Precompute the numerical flux for the exact equilibrium state ONCE
F_eq = compute_wb_flux(E_center)

def get_wb_afd_weno_rhs(u):
    # 1. Get numerical flux for the current state
    F_num = compute_wb_flux(u)
    
    # 2. Compute local scalar w
    w = u / E_center
    
    # 3. Discrete flux gradient
    dFdx = (F_num - np.roll(F_num, 1)) / dx
    
    # 4. Well-balanced discrete source term (S_i = w_i * S_{eq, i})
    S_discrete = w * (F_eq - np.roll(F_eq, 1)) / dx
    
    # 5. Final right-hand side
    dudt = -dFdx + S_discrete
    return dudt

# ---------------------------------------------------------
# 3. Main Simulation Loop
# ---------------------------------------------------------
# Set initial condition exactly to the steady-state
u = E_center.copy() 
u_init = u.copy()

# ---- VERIFICATION FOR THESIS ----
# If the scheme is well-balanced, dudt should be exactly zero.
# Note: Because the domain is periodic but the exponential is not, there is an 
# artificial shock at the boundary (index 0). We ignore cell 0 to check the interior.
initial_dudt = get_wb_afd_weno_rhs(u)
print(f"Max RHS error in interior at t=0: {np.max(np.abs(initial_dudt[1:])):.3e}")
# ---------------------------------

t = 0.0
while t < t_end:
    dt = CFL * dx / np.abs(c_speed)
    if t + dt > t_end:
        dt = t_end - t
        
    # SSP-RK3 Time Stepping
    u1 = u + dt * get_wb_afd_weno_rhs(u)
    u2 = 0.75 * u + 0.25 * u1 + 0.25 * dt * get_wb_afd_weno_rhs(u1)
    u = (1.0/3.0) * u + (2.0/3.0) * u2 + (2.0/3.0) * dt * get_wb_afd_weno_rhs(u2)
    
    t += dt

# ---------------------------------------------------------
# 4. Plot results
# ---------------------------------------------------------
plt.figure(figsize=(8, 5))
plt.plot(x, u_init, '--', label='Exact Equilibrium (t=0)', color='gray', linewidth=2)
plt.plot(x, u, '-', label=f'WB-AFD-WENO (t={t_end})', color='blue')
plt.title(f'Well-Balanced Advection-Reaction ($\\lambda$={lam})')
plt.xlabel('x')
plt.ylabel('u')
plt.legend()
plt.grid(True)
plt.show()
