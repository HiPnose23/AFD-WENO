#ifndef WENO_LIB_H_
#define WENO_LIB_H_

#include <cmath>
// =======================================================================
//  AMReX vs Vanilla C++ compatibility block 
// ======================================================================= 
#ifdef USE_AMREX
    #include <AMReX.H>
    #include <AMReX_REAL.H>
    using namespace amrex;
#else
    #define AMREX_GPU_DEVICE 
    #define AMREX_FORCE_INLINE inline
    using Real = double;
#endif

// ========================================================================
// WENO-AO(3) Interpolation
// ========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void weno_ao_3_interpolation(Real um2, Real um1, Real u0, Real up1, Real up2, 
                             Real dx, Real& u_L, Real& u_R, Real& dudx) 
{
    // Left-biased stencil S1
    Real pt1 = (25.0*u0 - 2.0*um1 + um2) / 24.0;
    Real x1_ = (3.0*u0 - 4.0*um1 + um2) / 2.0;
    Real x21 = (u0 - 2.0*um1 + um2) / 2.0;
    
    // Centered stencil S2
    Real pt2 = (22.0*u0 + um1 + up1) / 24.0;
    Real x2_ = (up1 - um1) / 2.0;
    Real x22 = (-2.0*u0 + um1 + up1) / 2.0;
    
    // Right-biased stencil S3
    Real pt3 = (25.0*u0 - 2.0*up1 + up2) / 24.0;
    Real x3_ = (-3.0*u0 + 4.0*up1 - up2) / 2.0;
    Real x23 = (u0 - 2.0*up1 + up2) / 2.0;
    
    // Smoothness indicators
    Real beta1 = x1_*x1_ + (13.0/3.0)*x21*x21;
    Real beta2 = x2_*x2_ + (13.0/3.0)*x22*x22;
    Real beta3 = x3_*x3_ + (13.0/3.0)*x23*x23;
    
    // Linear weights
    Real gamma_lo = 0.85;
    Real g1 = (1.0 - gamma_lo) / 2.0;
    Real g2 = gamma_lo;
    Real g3 = g1;
    
    // Non-linear weights mapping
    Real eps = 1e-12;
    Real tau = 0.5 * (std::abs(beta2 - beta1) + std::abs(beta2 - beta3));
    
    Real tmp1 = tau / (beta1 + eps);
    Real tmp2 = tau / (beta2 + eps);
    Real tmp3 = tau / (beta3 + eps);

    Real w1 = g1 * (1.0 + tmp1*tmp1);
    Real w2 = g2 * (1.0 + tmp2*tmp2);
    Real w3 = g3 * (1.0 + tmp3*tmp3);
    
    Real w_sum = w1 + w2 + w3;
    Real wb1 = w1 / w_sum;
    Real wb2 = w2 / w_sum;
    Real wb3 = w3 / w_sum;
    
    // Interpolate to left face (-1/2) and right face (+1/2)
    Real P1_R = pt1 + x1_*0.5 + x21*(1.0/6.0);
    Real P2_R = pt2 + x2_*0.5 + x22*(1.0/6.0);
    Real P3_R = pt3 + x3_*0.5 + x23*(1.0/6.0);
    
    Real P1_L = pt1 - x1_*0.5 + x21*(1.0/6.0);
    Real P2_L = pt2 - x2_*0.5 + x22*(1.0/6.0);
    Real P3_L = pt3 - x3_*0.5 + x23*(1.0/6.0);
    
    u_R = wb1*P1_R + wb2*P2_R + wb3*P3_R;
    u_L = wb1*P1_L + wb2*P2_L + wb3*P3_L;
    
    // Final derivative scaled by dx
    dudx = (wb1*x1_ + wb2*x2_ + wb3*x3_) / dx;
}

// ========================================================================
// WENO-AO(5,3) Interpolation
// ========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void weno_ao_5_3_interpolation(Real um2, Real um1, Real u0, Real up1, Real up2, 
                               Real dx, Real& u_L, Real& u_R, Real& dudx) 
{
    // r=3 stencils
    Real pt1=(25.0*u0-2.0*um1+um2)/24.0; Real x1_=(3.0*u0-4.0*um1+um2)/2.0; Real x21=(u0-2.0*um1+um2)/2.0;
    Real pt2=(22.0*u0+um1+up1)/24.0;     Real x2_=(up1-um1)/2.0;            Real x22=(-2.0*u0+um1+up1)/2.0;
    Real pt3=(25.0*u0-2.0*up1+up2)/24.0; Real x3_=(-3.0*u0+4.0*up1-up2)/2.0; Real x23=(u0-2.0*up1+up2)/2.0;
    
    Real b1 = x1_*x1_ + (13.0/3.0)*x21*x21; 
    Real b2 = x2_*x2_ + (13.0/3.0)*x22*x22; 
    Real b3 = x3_*x3_ + (13.0/3.0)*x23*x23;

    // r=5 central stencil
    Real upt5 = (5178.0*u0 + 308.0*(um1+up1) - 17.0*(um2+up2)) / 5760.0;
    Real ux5  = (-154.0*um1 + 17.0*um2 + 154.0*up1 - 17.0*up2) / 240.0;
    Real ux25 = (-402.0*u0 + 212.0*(um1+up1) - 11.0*(um2+up2)) / 336.0;
    Real ux35 = (2.0*um1 - um2 - 2.0*up1 + up2) / 12.0;
    Real ux45 = (6.0*u0 - 4.0*(um1+up1) + (um2+up2)) / 24.0;

    Real tmp5 = (ux25 + 123.0*ux45/455.0);
    Real b5 = ((ux5 + ux35/10.0)*(ux5 + ux35/10.0)
        + (13.0/3.0)*tmp5*tmp5
        + (781.0/20.0)*ux35*ux35
        + (1421461.0/2275.0)*ux45*ux45 / 1e10);

    // Linear weights
    Real gHi = 0.85; Real gLo = 0.85;
    Real g5  = gHi;
    Real g2  = (1.0 - gHi) * gLo;
    Real g1  = (1.0 - gHi) * (1.0 - gLo) / 2.0;
    Real g3  = g1;

    // tau
    Real eps  = 1e-12;
    Real tau5 = (std::abs(b5-b1) + std::abs(b5-b2) + std::abs(b5-b3)) / 3.0;
    Real tau3 = 0.5*(std::abs(b2-b1) + std::abs(b2-b3));

    // Nonlinear weights
    Real t5 = tau5/(b5+eps); Real w5 = g5 * (1.0 + t5*t5);
    Real t1 = tau3/(b1+eps); Real w1 = g1 * (1.0 + t1*t1);
    Real t2 = tau3/(b2+eps); Real w2 = g2 * (1.0 + t2*t2);
    Real t3 = tau3/(b3+eps); Real w3 = g3 * (1.0 + t3*t3);

    Real ws = w5+w1+w2+w3;
    Real wb5 = w5/ws; Real wb1 = w1/ws; Real wb2 = w2/ws; Real wb3 = w3/ws;

    // Face values
    Real P5R = upt5 + ux5*0.5  + ux25/6.0 + ux35/20.0  + ux45/70.0;
    Real P5L = upt5 - ux5*0.5  + ux25/6.0 - ux35/20.0  + ux45/70.0;
    Real P1R = pt1 + x1_*0.5 + x21/6.0; Real P2R = pt2 + x2_*0.5 + x22/6.0; Real P3R = pt3 + x3_*0.5 + x23/6.0;
    Real P1L = pt1 - x1_*0.5 + x21/6.0; Real P2L = pt2 - x2_*0.5 + x22/6.0; Real P3L = pt3 - x3_*0.5 + x23/6.0;

    u_R = (wb5/g5)*(P5R - g1*P1R - g2*P2R - g3*P3R) + wb1*P1R + wb2*P2R + wb3*P3R;
    u_L = (wb5/g5)*(P5L - g1*P1L - g2*P2L - g3*P3L) + wb1*P1L + wb2*P2L + wb3*P3L;

    dudx = ((wb5/g5)*((ux5 - (3.0/20.0)*ux35) - g1*x1_ - g2*x2_ - g3*x3_)
            + wb1*x1_ + wb2*x2_ + wb3*x3_) / dx;
}

// ========================================================================
// WENO-AO(4,3) Boundary Derivative
// ========================================================================
// Evaluated at face i+1/2 using f_centers from {i-1, i, i+1, i+2}
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void weno_ao_43_boundary(Real fm1, Real f0, Real f1, Real f2, Real dx, 
                         Real& d1dx, Real& d3dx) 
{
    Real eps = 1e-36;
    Real Hi = 0.85;

    // Large r=4 centered stencil Sr4zb
    Real fpt_c = (13.0*f0 - fm1 + 13.0*f1 - f2) / 24.0;
    Real fx_c  = (-63.0*f0 + fm1 + 63.0*f1 - f2) / 60.0;
    Real fx2_c = (-f0 + fm1 - f1 + f2) / 4.0;
    Real fx3_c = (3.0*f0 - fm1 - 3.0*f1 + f2) / 6.0;

    Real beta_c = (fx_c + fx3_c/10.0)*(fx_c + fx3_c/10.0) 
                  + (13.0/3.0)*fx2_c*fx2_c 
                  + (781.0/20.0)*fx3_c*fx3_c;

    // Small left-biased Sr3zb1
    Real fx1  = f1 - f0;
    Real fx21 = (f1 - 2.0*f0 + fm1) / 2.0;

    // Small right-biased Sr3zb2
    Real fx2  = f1 - f0;
    Real fx22 = (f2 - 2.0*f1 + f0) / 2.0;

    Real beta1 = fx1*fx1 + (13.0/3.0)*fx21*fx21;
    Real beta2 = fx2*fx2 + (13.0/3.0)*fx22*fx22;

    Real g_c = Hi;
    Real g1  = (1.0 - Hi) / 2.0;
    Real g2  = (1.0 - Hi) / 2.0;

    Real tau = (std::abs(beta_c - beta1) + std::abs(beta_c - beta2)) / 2.0;
    
    Real tc = tau/(beta_c + eps); Real wc = g_c * (1.0 + tc*tc);
    Real t1 = tau/(beta1  + eps); Real w1 = g1  * (1.0 + t1*t1);
    Real t2 = tau/(beta2  + eps); Real w2 = g2  * (1.0 + t2*t2);
    
    Real ws = wc + w1 + w2;
    wc /= ws; w1 /= ws; w2 /= ws;

    Real d1_c  = fx_c  + fx3_c * (-3.0/20.0);
    Real d1_1  = fx1;
    Real d1_2  = fx2;
    Real d3_c  = 6.0 * fx3_c;

    d1dx = ((wc/g_c)*(d1_c - g1*d1_1 - g2*d1_2) + w1*d1_1 + w2*d1_2) / dx;
    d3dx = ((wc/g_c) * d3_c) / dx;
}

#endif
