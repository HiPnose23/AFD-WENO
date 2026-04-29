# weno_lib.py
import numpy as np

def weno_ao_3_interpolation(u, dx, boundary_mode='wrap'):
    """
    Pointwise WENO-AO(3) interpolation.
    boundary_mode: 'wrap' for periodic, 'edge' for outflow/Sod/Lax.
    """
    # 1. Determine padding based on input shape (Scalar vs Vector)
    if u.ndim == 1:
        pad_width = (2, 2)
    else:
        pad_width = ((2, 2), (0, 0))
        
    U_pad = np.pad(u, pad_width, mode=boundary_mode)

        # Extract stencils relative to zone "0" (Eq 16)
    u_m2 = U_pad[0:-4]
    u_m1 = U_pad[1:-3]
    u_0  = U_pad[2:-2]
    u_1  = U_pad[3:-1]
    u_2  = U_pad[4:]
    
    # Left-biased stencil S1 (Eq 17)
    u_pt1 = (25.0*u_0 - 2.0*u_m1 + u_m2) / 24.0
    u_x1 = (3.0*u_0 - 4.0*u_m1 + u_m2) / 2.0
    u_x21 = (u_0 - 2.0*u_m1 + u_m2) / 2.0
    
    # Centered stencil S2 (Eq 18)
    u_pt2 = (22.0*u_0 + u_m1 + u_1) / 24.0
    u_x2 = (u_1 - u_m1) / 2.0
    u_x22 = (-2.0*u_0 + u_m1 + u_1) / 2.0
    
    # Right-biased stencil S3 (Eq 19)
    u_pt3 = (25.0*u_0 - 2.0*u_1 + u_2) / 24.0
    u_x3 = (-3.0*u_0 + 4.0*u_1 - u_2) / 2.0
    u_x23 = (u_0 - 2.0*u_1 + u_2) / 2.0
    
    # Smoothness indicators (Eq 20)
    beta1 = u_x1**2 + (13.0/3.0)*u_x21**2
    beta2 = u_x2**2 + (13.0/3.0)*u_x22**2
    beta3 = u_x3**2 + (13.0/3.0)*u_x23**2
    
    # Linear weights (Eq 21), using typical gamma_Lo = 0.85
    gamma_lo = 0.85
    g1 = (1.0 - gamma_lo) / 2.0
    g2 = gamma_lo
    g3 = (1.0 - gamma_lo) / 2.0
    
    # Non-linear weights mapping (Eq 22, 23, 24)
    epsilon = 1e-12
    tau = 0.5 * (np.abs(beta2 - beta1) + np.abs(beta2 - beta3))
    
    w1 = g1 * (1.0 + (tau / (beta1 + epsilon))**2)
    w2 = g2 * (1.0 + (tau / (beta2 + epsilon))**2)
    w3 = g3 * (1.0 + (tau / (beta3 + epsilon))**2)
    
    w_sum = w1 + w2 + w3
    wb1 = w1 / w_sum
    wb2 = w2 / w_sum
    wb3 = w3 / w_sum
    
    # Interpolate to left face (-1/2) and right face (+1/2) using Legendre polynomials
    # P(x) = u_pt + u_x * L1(x) + u_x2 * L2(x)
    # At x = 1/2: L1=1/2, L2=1/6
    # At x =-1/2: L1=-1/2, L2=1/6
    P1_R = u_pt1 + u_x1*0.5 + u_x21*(1.0/6.0)
    P2_R = u_pt2 + u_x2*0.5 + u_x22*(1.0/6.0)
    P3_R = u_pt3 + u_x3*0.5 + u_x23*(1.0/6.0)
    
    P1_L = u_pt1 - u_x1*0.5 + u_x21*(1.0/6.0)
    P2_L = u_pt2 - u_x2*0.5 + u_x22*(1.0/6.0)
    P3_L = u_pt3 - u_x3*0.5 + u_x23*(1.0/6.0)
    
    U_R_eval = wb1*P1_R + wb2*P2_R + wb3*P3_R
    U_L_eval = wb1*P1_L + wb2*P2_L + wb3*P3_L
    
    # Final derivative scaled by dx
    dUdx = (wb1*u_x1 + wb2*u_x2 + wb3*u_x3) / dx
    
    return U_L_eval, U_R_eval, dUdx

def weno_ao_5_3_interpolation(u, dx, boundary_mode='wrap'):
    """
    WENO-AO(5,3): non-linear hybridization of one large r=5 central stencil
    with three small r=3 CWENO stencils (Eqs. 26–33, Balsara et al. 2023).
    Falls back to stable 3rd-order CWENO near shocks/discontinuities.
    """
    pad_width = (2, 2) if u.ndim == 1 else ((2, 2), (0, 0))
    U_pad = np.pad(u, pad_width, mode=boundary_mode)
    u_m2 = U_pad[0:-4]; u_m1 = U_pad[1:-3]
    u_0  = U_pad[2:-2]; u_1  = U_pad[3:-1]; u_2 = U_pad[4:]

    # === r=3 stencils (Eqs. 17–20) ===
    pt1=(25*u_0-2*u_m1+u_m2)/24; x1_=(3*u_0-4*u_m1+u_m2)/2; x21=(u_0-2*u_m1+u_m2)/2
    pt2=(22*u_0+u_m1+u_1)/24;    x2_=(u_1-u_m1)/2;           x22=(-2*u_0+u_m1+u_1)/2
    pt3=(25*u_0-2*u_1+u_2)/24;   x3_=(-3*u_0+4*u_1-u_2)/2;  x23=(u_0-2*u_1+u_2)/2
    b1=x1_**2+(13/3)*x21**2; b2=x2_**2+(13/3)*x22**2; b3=x3_**2+(13/3)*x23**2

    # === r=5 central stencil (Eq. 27) ===
    upt5 = (5178*u_0 + 308*(u_m1+u_1) - 17*(u_m2+u_2)) / 5760
    ux5  = (-154*u_m1 + 17*u_m2 + 154*u_1 - 17*u_2) / 240
    ux25 = (-402*u_0 + 212*(u_m1+u_1) - 11*(u_m2+u_2)) / 336
    ux35 = (2*u_m1 - u_m2 - 2*u_1 + u_2) / 12
    ux45 = (6*u_0 - 4*(u_m1+u_1) + (u_m2+u_2)) / 24

    # Eq. 28: smoothness indicator for r=5 stencil
    b5 = ((ux5 + ux35/10)
        + (13/3)*(ux25 + 123*ux45/455)**2
        + (781/20)*ux35**2
        + (1421461/2275)*ux45**2 / (1e10)) 

    # === Linear weights (Eq. 29) ===
    gHi = 0.85; gLo = 0.85
    g5  = gHi
    g2  = (1 - gHi) * gLo
    g1  = (1 - gHi) * (1 - gLo) / 2
    g3  = (1 - gHi) * (1 - gLo) / 2

    # === τ for the large stencil vs all three small stencils (Eq. 30) ===
    eps  = 1e-12
    tau5 = (np.abs(b5-b1) + np.abs(b5-b2) + np.abs(b5-b3)) / 3
    tau3 = 0.5*(np.abs(b2-b1) + np.abs(b2-b3))  # used for r=3 sub-weights

    # === Unnormalized nonlinear weights (Eq. 31) ===
    w5 = g5 * (1 + (tau5/(b5+eps))**2)
    w1 = g1 * (1 + (tau3/(b1+eps))**2)
    w2 = g2 * (1 + (tau3/(b2+eps))**2)
    w3 = g3 * (1 + (tau3/(b3+eps))**2)

    # === Normalization (Eq. 32) ===
    ws = w5+w1+w2+w3
    wb5=w5/ws; wb1=w1/ws; wb2=w2/ws; wb3=w3/ws

    # === Face values via Legendre basis ===
    # L1(±½)=±½, L2(±½)=1/6, L3(±½)=±1/20, L4(±½)=1/70
    P5R = upt5 + ux5*0.5  + ux25/6 + ux35/20  + ux45/70
    P5L = upt5 - ux5*0.5  + ux25/6 - ux35/20  + ux45/70
    P1R=pt1+x1_*0.5+x21/6; P2R=pt2+x2_*0.5+x22/6; P3R=pt3+x3_*0.5+x23/6
    P1L=pt1-x1_*0.5+x21/6; P2L=pt2-x2_*0.5+x22/6; P3L=pt3-x3_*0.5+x23/6

    # === WENO-AO(5,3) combination (Eq. 33) ===
    # In smooth limit: wb5→g5, so (wb5/g5)*(P5 - g1*P1-g2*P2-g3*P3)
    #                            + wb1*P1+wb2*P2+wb3*P3 → P5  ✓
    UR = (wb5/g5)*(P5R - g1*P1R - g2*P2R - g3*P3R) + wb1*P1R + wb2*P2R + wb3*P3R
    UL = (wb5/g5)*(P5L - g1*P1L - g2*P2L - g3*P3L) + wb1*P1L + wb2*P2L + wb3*P3L

    # Derivative at cell center (ux coefficient from each stencil)
    dUdx = ((wb5/g5)*(ux5 - g1*x1_ - g2*x2_ - g3*x3_)
            + wb1*x1_ + wb2*x2_ + wb3*x3_) / dx

    return UL, UR, dUdx

def weno_ao_43_boundary(f, dx, Hi=0.85, eps=1e-36, boundary_mode = 'wrap'):
    """
    Section 4.2. Input f[i] = (A * ux)[i] at zone centers.
    Output at each face i+1/2:
      d1[i] = dP/dx|_{x=0}  = fx  - (3/20)*fx3       (physical: divide by dx)
      d3[i] = d³P/dx³|_{x=0} = 6 * fx3               (physical: divide by dx³)
    """
    pad = (1, 2) if f.ndim == 1 else ((1, 2), (0, 0))
    f_pad = np.pad(f, pad, mode=boundary_mode)
    fm1 = f_pad[0:-3]
    f0  = f_pad[1:-2]
    f1  = f_pad[2:-1]
    f2  = f_pad[3:]
    
    # === Large r=4 centered stencil Sr4zb (Eq. 54) ===
    # Origin at face i+1/2, between f0 and f1
    fpt_c = (13*f0 - fm1 + 13*f1 - f2) / 24
    fx_c  = (-63*f0 + fm1 + 63*f1 - f2) / 60  # L1 coeff
    fx2_c = (-f0 + fm1 - f1 + f2) / 4          # L2 coeff
    fx3_c = (3*f0 - fm1 - 3*f1 + f2) / 6       # L3 coeff

    # Smoothness indicator for r=4c (from BGS16 structure)
    beta_c = (fx_c + fx3_c/10)**2 + (13/3)*(fx2_c)**2 + (781/20)*fx3_c**2

    # === Small left-biased Sr3zb1: uses fm1, f0, f1 (Eq. 48) ===
    fpt1 = (8*f0 - fm1 + 5*f1) / 12
    fx1  = f1 - f0
    fx21 = (f1 - 2*f0 + fm1) / 2

    # === Small right-biased Sr3zb2: uses f0, f1, f2 (Eq. 49) ===
    fpt2 = (8*f1 - f2 + 5*f0) / 12
    fx2  = f1 - f0
    fx22 = (f2 - 2*f1 + f0) / 2

    beta1 = fx1**2 + (13/3)*fx21**2
    beta2 = fx2**2 + (13/3)*fx22**2

    # === Linear weights (Eq. 55) ===
    g_c = Hi
    g1  = (1 - Hi) / 2
    g2  = (1 - Hi) / 2

    tau = (abs(beta_c - beta1) + abs(beta_c - beta2)) / 2
    wc = g_c * (1 + (tau/(beta_c + eps))**2)
    w1 = g1  * (1 + (tau/(beta1  + eps))**2)
    w2 = g2  * (1 + (tau/(beta2  + eps))**2)
    ws = wc + w1 + w2
    wc /= ws; w1 /= ws; w2 /= ws

    # === AO combination (Eq. 57) ===
    # PAO = (wc/gc)*(Pc - g1*P1 - g2*P2) + w1*P1 + w2*P2

    # Evaluate derivatives at x=0 (the face):
    # dL1/dx|_0 = 1,  dL2/dx|_0 = 0,  dL3/dx|_0 = -3/20
    # d³L3/dx³  = 6,  all others zero at any x
    d1_c  = fx_c  + fx3_c * (-3/20)     # dPc/dx at x=0
    d1_1  = fx1   + 0                   # dP1/dx at x=0 (only L1 and L2 terms)
    d1_2  = fx2   + 0

    d3_c  = 6 * fx3_c                   # d³Pc/dx³ at x=0
    # small stencils only go to L2, so d³ = 0

    d1dx = ((wc/g_c)*(d1_c - g1*d1_1 - g2*d1_2) + w1*d1_1 + w2*d1_2) / dx
    d3dx = ((wc/g_c) * d3_c) / dx   # small stencil contributions cancel out

    return d1dx, d3dx
