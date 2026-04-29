import numpy as np
import matplotlib.pyplot as plt
from weno_lib import *

# =========================================================
# 1. Physics Setup (Euler with Gravity)
# =========================================================
gamma = 1.4
nx = 200
x, dx = np.linspace(-1.0, 1.0, nx, endpoint=False, retstep=True)

# Gravity potential phi(x) = x (meaning g = 1 pulling left)
# Isothermal Equilibrium state: p_eq = rho_eq = exp(-phi)
rho_eq_c = np.exp(-x)
p_eq_c   = np.exp(-x)

x_face = x + dx / 2.0
rho_eq_f = np.exp(-x_face)
p_eq_f   = np.exp(-x_face)

# Analytical derivatives of the equilibrium state
rho_eq_prime_c = -rho_eq_c
p_eq_prime_c   = -p_eq_c

# =========================================================
# 2. Variable Transformations
# =========================================================
def U_to_W(U, rho_e, p_e):
    """Maps Conservative Variables U to Equilibrium Variables W"""
    rho = U[0]
    rhou = U[1]
    E = U[2]
    
    u = rhou / rho
    p = (gamma - 1.0) * (E - 0.5 * rho * u**2)
    
    w1 = rho / rho_e
    w2 = u
    w3 = p / p_e
    return np.array([w1, w2, w3]), p

def W_to_U(W, rho_e, p_e):
    """Maps Equilibrium Variables W back to Conservative Variables U"""
    rho = W[0] * rho_e
    u = W[1]
    p = W[2] * p_e
    
    rhou = rho * u
    E = p / (gamma - 1.0) + 0.5 * rho * u**2
    return np.array([rho, rhou, E])

def euler_flux(U, p):
    """Computes the standard physical Euler flux vector"""
    rho, rhou, E = U[0], U[1], U[2]
    u = rhou / rho
    
    F1 = rhou
    F2 = rhou * u + p
    F3 = u * (E + p)
    return np.array([F1, F2, F3])

# =========================================================
# 3. Well-Balanced AFD-WENO Fluxes
# =========================================================
def compute_wb_flux(U):
    """Computes 3rd-order WB AFD-WENO fluxes for Euler equations"""
    W, p_c = U_to_W(U, rho_eq_c, p_eq_c)
    
    W_L = np.zeros_like(W)
    W_R = np.zeros_like(W)
    Wx_c = np.zeros_like(W)
    
    # 1. WENO-AO(3) on Equilibrium Variables
    for k in range(3):
        w_L, w_R, dwdx = weno_ao_3_interpolation(W[k], dx, 'wrap')
        W_L[k] = w_R
        W_R[k] = np.roll(w_L, -1)
        Wx_c[k] = dwdx
        
    # 2. Reconstruct U at faces
    U_L_face = W_to_U(W_L, rho_eq_f, p_eq_f)
    U_R_face = W_to_U(W_R, rho_eq_f, p_eq_f)
    
    # 3. Local Lax-Friedrichs (Rusanov) Riemann Solver
    _, p_L = U_to_W(U_L_face, rho_eq_f, p_eq_f)
    _, p_R = U_to_W(U_R_face, rho_eq_f, p_eq_f)
    
    c_L = np.sqrt(gamma * p_L / U_L_face[0])
    c_R = np.sqrt(gamma * p_R / U_R_face[0])
    u_L = U_L_face[1] / U_L_face[0]
    u_R = U_R_face[1] / U_R_face[0]
    
    # Max wave speed at each face (broadcast to shape 1, nx for matrix math)
    alpha = np.maximum(np.abs(u_L) + c_L, np.abs(u_R) + c_R).reshape(1, nx)
    
    F_L = euler_flux(U_L_face, p_L)
    F_R = euler_flux(U_R_face, p_R)
    F_star = 0.5 * (F_L + F_R - alpha * (U_R_face - U_L_face))
    
    # 4. Center derivatives (Chain Rule mapping Wx back to Ux)
    w1, w2, w3 = W[0], W[1], W[2]
    w1x, w2x, w3x = Wx_c[0], Wx_c[1], Wx_c[2]
    
    rho_x = w1x * rho_eq_c + w1 * rho_eq_prime_c
    u_x   = w2x
    p_x   = w3x * p_eq_c   + w3 * p_eq_prime_c
    
    rho, u, p, E = U[0], U[1]/U[0], p_c, U[2]
    
    # Differentiating the physical flux F(U)
    rhou_x = rho_x * u + rho * u_x
    rhou2_p_x = u * rhou_x + rho * u * u_x + p_x
    uEp_x = u_x * (E + p) + u * (p_x / (gamma - 1.0) + 0.5 * (rho_x * u**2 + 2*rho*u*u_x) + p_x)
    
    fc = np.array([rhou_x, rhou2_p_x, uEp_x])
    
    # 5. AFD-WENO flux correction (Strictly 3rd order)
    F_num = np.zeros_like(F_star)
    for k in range(3):
        d1dx, _ = weno_ao_43_boundary(fc[k], dx, boundary_mode='wrap')
        F_num[k] = F_star[k] - (dx**2 / 24.0) * d1dx
        
    return F_num

# Precompute the numerical flux for the exact equilibrium state ONCE
# U_eq = [rho_eq, 0, E_eq]
U_eq = W_to_U(np.array([np.ones(nx), np.zeros(nx), np.ones(nx)]), rho_eq_c, p_eq_c)
F_eq = compute_wb_flux(U_eq)

def get_euler_wb_rhs(U):
    """Evaluates the full RHS with the Xing-Shu source term balance"""
    F_num = compute_wb_flux(U)
    
    # Flux divergence
    dFdx = (F_num - np.roll(F_num, 1, axis=1)) / dx
    
    # Well-balanced discrete source term
    W, _ = U_to_W(U, rho_eq_c, p_eq_c)
    w1, w2 = W[0], W[1]
    
    S_discrete = np.zeros_like(U)
    
    # Momentum source relies on the discrete divergence of equilibrium pressure flux
    F_eq_2_div = (F_eq[1] - np.roll(F_eq[1], 1)) / dx
    
    S_discrete[0] = 0.0
    S_discrete[1] = w1 * F_eq_2_div      # w1 * d(p_eq)/dx
    S_discrete[2] = w2 * S_discrete[1]   # u * Momentum Source
    
    return -dFdx + S_discrete

# =========================================================
# 4. Rigorous Testing Suite
# =========================================================
def run_tests():
    print("\n--- TEST 1: EXACT STEADY-STATE PRESERVATION ---")
    # Initialize with exact isothermal equilibrium
    U = U_eq.copy()
    rhs = get_euler_wb_rhs(U)
    
    interior = slice(5, -5) # Ignore boundaries polluted by periodic wrapping of exponential
    max_err_rho = np.max(np.abs(rhs[0, interior]))
    max_err_mom = np.max(np.abs(rhs[1, interior]))
    max_err_en  = np.max(np.abs(rhs[2, interior]))
    
    print(f"Density RHS Error:  {max_err_rho:.3e}")
    print(f"Momentum RHS Error: {max_err_mom:.3e}")
    print(f"Energy RHS Error:   {max_err_en:.3e}")
    
    print("\n--- TEST 2: SMALL PERTURBATION ADVECTION ---")
    print("Simulating acoustic wave on top of background gravity...")
    # Add a tiny pressure perturbation (1e-5) to the equilibrium state
    W_init, _ = U_to_W(U_eq, rho_eq_c, p_eq_c)
    W_init[2] += 1e-5 * np.exp(-200 * (x + 0.2)**2) # Tiny pressure pulse
    U_pert = W_to_U(W_init, rho_eq_c, p_eq_c)
    
    U = U_pert.copy()
    t, t_end = 0.0, 0.15 # Run short time to keep periodic wrap shock away from center
    CFL = 0.4
    
    while t < t_end:
        # Compute max timestep based on acoustic speed
        _, p = U_to_W(U, rho_eq_c, p_eq_c)
        c = np.sqrt(gamma * p / U[0])
        max_speed = np.max(np.abs(U[1]/U[0]) + c)
        dt = CFL * dx / max_speed
        if t + dt > t_end: dt = t_end - t
        
        # SSP-RK3
        U1 = U + dt * get_euler_wb_rhs(U)
        U2 = 0.75 * U + 0.25 * U1 + 0.25 * dt * get_euler_wb_rhs(U1)
        U  = (1/3) * U + (2/3) * U2 + (2/3) * dt * get_euler_wb_rhs(U2)
        t += dt

    # Plot the Pressure Perturbation (Current Pressure - Equilibrium Pressure)
    _, p_final = U_to_W(U, rho_eq_c, p_eq_c)
    p_pert_final = p_final - p_eq_c
    
    plt.figure(figsize=(9, 5))
    plt.plot(x, p_pert_final, 'b-', linewidth=2, label=f'WB-AFD-WENO (t={t_end})')
    plt.title('Evolution of $10^{-5}$ Pressure Perturbation in Gravitational Field')
    plt.xlabel('x')
    plt.ylabel('Pressure Perturbation $\\Delta p$')
    plt.xlim(-0.5, 0.5)
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    run_tests()
