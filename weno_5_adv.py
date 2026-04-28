import numpy as np
import matplotlib.pyplot as plt
from weno_lib import *
# ---------------------------------------------------------
# 1. Setup and Parameters (Advection-Reaction)
# ---------------------------------------------------------
c_speed = 1.0    # Advection speed (c)
lam = -0.2       # Reaction rate (lambda) -> Negative is decay
CFL = 0.8
t_end = 2.0      # End time (one full period if c=1 and L=2)
nx = 200
x = np.linspace(-1.0, 1.0, nx, endpoint=False) # Periodic domain [-1, 1]
dx = x[1] - x[0]


# ---------------------------------------------------------
# 2/3. AFD-WENO RHS (Advection-Reaction Version)
# ---------------------------------------------------------
def get_afd_weno_rhs(u):
    u_L_eval, u_R_eval, dudx_center = weno_ao_3_interpolation(u, dx, 'wrap')
    
    # Boundary states at faces i+1/2 (with periodic wrap)
    u_L_face = u_R_eval
    u_R_face = np.roll(u_L_eval, -1)
    
    # Step (v): Upwind Riemann Flux (since c is constant)
    if c_speed >= 0:
        F_star = c_speed * u_L_face
    else:
        F_star = c_speed * u_R_face
        
    # Step (vi): f = A * dudx. For advection, A is just the scalar 'c'
    f_centers = c_speed * dudx_center
    
    # Step (vii): Boundary derivative of f
    # For periodic, we roll f to get the difference across the boundary
    d1dx, d3dx = weno_ao_43_boundary(f_centers, dx)
    # High order AFD-WENO flux correction (Eq 14)
    F_num = F_star - (dx**2 / 24.0) * d1dx + (7 * dx**4 / 5760) * d3dx
    
    # Final RHS: -dF/dx + lambda * u
    # Note: F_num[i] is the flux at face i+1/2
    dudt = -(F_num - np.roll(F_num, 1)) / dx + lam * u
    return dudt

# ---------------------------------------------------------
# 4. Main Simulation Loop (Smooth Sine Wave)
# ---------------------------------------------------------
# Initial condition: Sine wave
u = np.sin(np.pi * x)
u_init = u.copy() # Store for plotting

t = 0.0
while t < t_end:
    dt = CFL * dx / np.abs(c_speed)
    if t + dt > t_end:
        dt = t_end - t
        
    # SSP-RK3 Time Stepping
    u1 = u + dt * get_afd_weno_rhs(u)
    u2 = 0.75 * u + 0.25 * u1 + 0.25 * dt * get_afd_weno_rhs(u1)
    u = (1.0/3.0) * u + (2.0/3.0) * u2 + (2.0/3.0) * dt * get_afd_weno_rhs(u2)
    
    t += dt

# Plot results
plt.figure(figsize=(8, 5))
plt.plot(x, u_init, '--', label='Initial (t=0)', color='gray')
plt.plot(x, u, 'o-', markersize=3, label=f'AFD-WENO (t={t_end})', color='blue')
plt.title(f'Advection-Reaction ($\lambda$={lam})')
plt.xlabel('x')
plt.ylabel('u')
plt.legend()
plt.grid(True)
plt.show()
