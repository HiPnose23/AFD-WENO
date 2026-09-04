#ifndef WENO_LIB_H_
#define WENO_LIB_H_

#include <cmath>


// ========================================================================
// WENO-AO(3) Interpolation
// ========================================================================
void weno_ao_3_interpolation(double um2, double um1, double u0, double up1, double up2, 
                             double dx, double& u_L, double& u_R, double& dudx) 
{
    // Left-biased stencil S1
    double pt1 = (25.0*u0 - 2.0*um1 + um2) / 24.0;
    double x1_ = (3.0*u0 - 4.0*um1 + um2) / 2.0;
    double x21 = (u0 - 2.0*um1 + um2) / 2.0;
    
    // Centered stencil S2
    double pt2 = (22.0*u0 + um1 + up1) / 24.0;
    double x2_ = (up1 - um1) / 2.0;
    double x22 = (-2.0*u0 + um1 + up1) / 2.0;
    
    // Right-biased stencil S3
    double pt3 = (25.0*u0 - 2.0*up1 + up2) / 24.0;
    double x3_ = (-3.0*u0 + 4.0*up1 - up2) / 2.0;
    double x23 = (u0 - 2.0*up1 + up2) / 2.0;
    
    // Smoothness indicators
    double beta1 = x1_*x1_ + (13.0/3.0)*x21*x21;
    double beta2 = x2_*x2_ + (13.0/3.0)*x22*x22;
    double beta3 = x3_*x3_ + (13.0/3.0)*x23*x23;
    
    // Linear weights
    double gamma_lo = 0.85;
    double g1 = (1.0 - gamma_lo) / 2.0;
    double g2 = gamma_lo;
    double g3 = g1;
    
    // Non-linear weights mapping
    double eps = 1e-12;
    double tau = 0.5 * (std::abs(beta2 - beta1) + std::abs(beta2 - beta3));
    
    double tmp1 = tau / (beta1 + eps);
    double tmp2 = tau / (beta2 + eps);
    double tmp3 = tau / (beta3 + eps);

    double w1 = g1 * (1.0 + tmp1*tmp1);
    double w2 = g2 * (1.0 + tmp2*tmp2);
    double w3 = g3 * (1.0 + tmp3*tmp3);
    
    double w_sum = w1 + w2 + w3;
    double wb1 = w1 / w_sum;
    double wb2 = w2 / w_sum;
    double wb3 = w3 / w_sum;
    
    // Interpolate to left face (-1/2) and right face (+1/2)
    double P1_R = pt1 + x1_*0.5 + x21*(1.0/6.0);
    double P2_R = pt2 + x2_*0.5 + x22*(1.0/6.0);
    double P3_R = pt3 + x3_*0.5 + x23*(1.0/6.0);
    
    double P1_L = pt1 - x1_*0.5 + x21*(1.0/6.0);
    double P2_L = pt2 - x2_*0.5 + x22*(1.0/6.0);
    double P3_L = pt3 - x3_*0.5 + x23*(1.0/6.0);
    
    u_R = wb1*P1_R + wb2*P2_R + wb3*P3_R;
    u_L = wb1*P1_L + wb2*P2_L + wb3*P3_L;
    
    // Final derivative scaled by dx
    dudx = (wb1*x1_ + wb2*x2_ + wb3*x3_) / dx;
}

// ========================================================================
// WENO-AO(5,3) Interpolation
// ========================================================================
void weno_ao_5_3_interpolation(double um2, double um1, double u0, double up1, double up2, 
                               double dx, double& u_L, double& u_R, double& dudx) 
{
    // r=3 stencils
    double pt1=(25.0*u0-2.0*um1+um2)/24.0; double x1_=(3.0*u0-4.0*um1+um2)/2.0; double x21=(u0-2.0*um1+um2)/2.0;
    double pt2=(22.0*u0+um1+up1)/24.0;     double x2_=(up1-um1)/2.0;            double x22=(-2.0*u0+um1+up1)/2.0;
    double pt3=(25.0*u0-2.0*up1+up2)/24.0; double x3_=(-3.0*u0+4.0*up1-up2)/2.0; double x23=(u0-2.0*up1+up2)/2.0;
    
    double b1 = x1_*x1_ + (13.0/3.0)*x21*x21; 
    double b2 = x2_*x2_ + (13.0/3.0)*x22*x22; 
    double b3 = x3_*x3_ + (13.0/3.0)*x23*x23;

    // r=5 central stencil
    double upt5 = (5178.0*u0 + 308.0*(um1+up1) - 17.0*(um2+up2)) / 5760.0;
    double ux5  = (-154.0*um1 + 17.0*um2 + 154.0*up1 - 17.0*up2) / 240.0;
    double ux25 = (-402.0*u0 + 212.0*(um1+up1) - 11.0*(um2+up2)) / 336.0;
    double ux35 = (2.0*um1 - um2 - 2.0*up1 + up2) / 12.0;
    double ux45 = (6.0*u0 - 4.0*(um1+up1) + (um2+up2)) / 24.0;

    double tmp5 = (ux25 + 123.0*ux45/455.0);
    double b5 = ((ux5 + ux35/10.0)*(ux5 + ux35/10.0)
        + (13.0/3.0)*tmp5*tmp5
        + (781.0/20.0)*ux35*ux35
        + (1421461.0/2275.0)*ux45*ux45);

    // Linear weights
    double gHi = 0.85; double gLo = 0.85; 
    double g5  = gHi;
    double g2  = (1.0 - gHi) * gLo;
    double g1  = (1.0 - gHi) * (1.0 - gLo) / 2.0;
    double g3  = g1;

    // tau
    double eps  = 1e-36;
    double tau5 = (std::abs(b5-b1) + std::abs(b5-b2) + std::abs(b5-b3)) / 3.0;
    double tau3 = 0.5*(std::abs(b2-b1) + std::abs(b2-b3));

    // Nonlinear weights
    double t5 = tau5/(b5+eps); double w5 = g5 * (1.0 + t5*t5);
    double t1 = tau5/(b1+eps); double w1 = g1 * (1.0 + t1*t1);
    double t2 = tau5/(b2+eps); double w2 = g2 * (1.0 + t2*t2);
    double t3 = tau5/(b3+eps); double w3 = g3 * (1.0 + t3*t3);

    double ws = w5+w1+w2+w3;
    double wb5 = w5/ws; double wb1 = w1/ws; double wb2 = w2/ws; double wb3 = w3/ws;

    // Face values
    double P5R = upt5 + ux5*0.5  + ux25/6.0 + ux35/20.0  + ux45/70.0;
    double P5L = upt5 - ux5*0.5  + ux25/6.0 - ux35/20.0  + ux45/70.0;
    double P1R = pt1 + x1_*0.5 + x21/6.0; double P2R = pt2 + x2_*0.5 + x22/6.0; double P3R = pt3 + x3_*0.5 + x23/6.0;
    double P1L = pt1 - x1_*0.5 + x21/6.0; double P2L = pt2 - x2_*0.5 + x22/6.0; double P3L = pt3 - x3_*0.5 + x23/6.0;

    u_R = (wb5/g5)*(P5R - g1*P1R - g2*P2R - g3*P3R) + wb1*P1R + wb2*P2R + wb3*P3R;
    u_L = (wb5/g5)*(P5L - g1*P1L - g2*P2L - g3*P3L) + wb1*P1L + wb2*P2L + wb3*P3L;

    dudx = ((wb5/g5)*((ux5 - (3.0/20.0)*ux35) - g1*x1_ - g2*x2_ - g3*x3_)
            + wb1*x1_ + wb2*x2_ + wb3*x3_) / dx;
}

// ========================================================================
// WENO-AO(4,3) Boundary Derivative
// ========================================================================
// Evaluated at face i+1/2 using f_centers from {i-1, i, i+1, i+2}
void weno_ao_43_boundary(double fm1, double f0, double f1, double f2, double dx, 
                         double& d1dx, double& d3dx, double& d2dx2) 
{
    double eps = 1e-36;
    double Hi = 0.85;

    // Large r=4 centered stencil Sr4zb
    double fpt_c = (13.0*f0 - fm1 + 13.0*f1 - f2) / 24.0;
    double fx_c  = (-63.0*f0 + fm1 + 63.0*f1 - f2) / 60.0;
    double fx2_c = (-f0 + fm1 - f1 + f2) / 4.0;
    double fx3_c = (3.0*f0 - fm1 - 3.0*f1 + f2) / 6.0;

    double beta_c = (fx_c + fx3_c/10.0)*(fx_c + fx3_c/10.0) 
                  + (13.0/3.0)*fx2_c*fx2_c 
                  + (781.0/20.0)*fx3_c*fx3_c;

    // Small left-biased Sr3zb1
    double fx1  = f1 - f0;
    double fx21 = (f1 - 2.0*f0 + fm1) / 2.0;

    // Small right-biased Sr3zb2
    double fx2  = f1 - f0;
    double fx22 = (f2 - 2.0*f1 + f0) / 2.0;

    double beta1 = fx1*fx1 + (13.0/3.0)*fx21*fx21;
    double beta2 = fx2*fx2 + (13.0/3.0)*fx22*fx22;

    double g_c = Hi;
    double g1  = (1.0 - Hi) / 2.0;
    double g2  = (1.0 - Hi) / 2.0;

    double tau = (std::abs(beta_c - beta1) + std::abs(beta_c - beta2)) / 2.0;
    
    double tc = tau/(beta_c + eps);
    double wc = g_c * (1.0 + tc*tc);
    double t1 = tau/(beta1  + eps);
    double w1 = g1  * (1.0 + t1*t1);
    double t2 = tau/(beta2  + eps);
    double w2 = g2  * (1.0 + t2*t2);
    
    double ws = wc + w1 + w2;
    wc /= ws; w1 /= ws; w2 /= ws;

    double d2_c  = 2.0 * fx2_c;
    double d2_1  = 2.0 * fx21;
    double d2_2  = 2.0 * fx22;

    d2dx2 = ((wc/g_c)*(d2_c - g1*d2_1 - g2*d2_2) + w1*d2_1 + w2*d2_2) / (dx*dx);
}

void weno_ao_63_boundary(double fm2, double fm1, double f0, double f1, double f2, double f3, 
                         double dx, double& d2dx2, double& d4dx4) 
{
    double eps = 1e-36;
    double Hi = 0.85; 

    // Variables to avoid writing coefficients twice 
    double sum0 = f1 + f0;
    double sum1 = f2 + fm1;
    double sum2 = f3 + fm2;

    double dif0 = f1 - f0;
    double dif1 = f2 - fm1;
    double dif2 = f3 - fm2;

    // Large r=6 centered stencil Sr6zbc (Equation 59)
    double fpt_c = (802.0*sum0 - 93.0*sum1 + 11.0*sum2) / 1440.0;
    double fx_c  = (1794.0*dif0 - 43.0*dif1 + 3.0*dif2) / 1680.0;
    double fx2_c = (-29.0*sum0 + 33.0*sum1 - 4.0*sum2) / 84.0;
    double fx3_c = (37.0*dif0 - 14.0*dif1 + dif2) / 54.0;
    double fx4_c = (2.0*sum0 - 3.0*sum1 + sum2) / 48.0;
    double fx5_c = (10.0*dif0 - 5.0*dif1 + dif2) / 120.0;
    
    double beta_c = (fx_c + fx3_c/10.0 + fx5_c/126.0)*(fx_c + fx3_c/10.0 + fx5_c/126.0)
                    + 13.0/3.0*(fx2_c+123.0/455.0*fx4_c)*(fx2_c+123.0/455.0*fx4_c)
                    + 781.0/20.0*(fx3_c+26045.0/49203.0*fx5_c)*(fx3_c+26045.0/49203.0*fx5_c)
                    + 1421461.0/2275.0*fx4_c*fx4_c
                    + 21520059541.0/1377684.0*fx5_c*fx5_c; 

    // Small left-biased Sr3zb1
    double fx1  = f1 - f0;
    double fx21 = (f1 - 2.0*f0 + fm1) / 2.0;

    // Small right-biased Sr3zb2
    double fx2  = f1 - f0;
    double fx22 = (f2 - 2.0*f1 + f0) / 2.0;

    double beta1 = fx1*fx1 + (13.0/3.0)*fx21*fx21;
    double beta2 = fx2*fx2 + (13.0/3.0)*fx22*fx22;

    double g_c = Hi;
    double g1  = (1.0 - Hi) / 2.0;
    double g2  = (1.0 - Hi) / 2.0;

    double tau = (std::abs(beta_c - beta1) + std::abs(beta_c - beta2)) / 2.0;
    
    double tc = tau / (beta_c + eps);
    double wc = g_c * (1.0 + tc*tc);
    double t1 = tau / (beta1  + eps);
    double w1 = g1  * (1.0 + t1*t1);
    double t2 = tau / (beta2  + eps);
    double w2 = g2  * (1.0 + t2*t2);
    
    double ws = wc + w1 + w2;
    wc /= ws; w1 /= ws; w2 /= ws;

    // --- EVEN DERIVATIVES FOR FLUX CORRECTIONS ---
    
    // The exact 2nd derivative of the 6th order polynomial Pr6zbc(x) evaluated at x=0
    // d^2/dx^2 [ fpt + fx*L1 + fx2*L2 + fx3*L3 + fx4*L4 + fx5*L5 ]
    // L2(0)'' = 2, L4(0)'' = -3/7
    double d2_c = 2.0 * fx2_c - (3.0 / 7.0) * fx4_c;
    double d2_1 = 2.0 * fx21;
    double d2_2 = 2.0 * fx22;

    // The exact 4th derivative of the 6th order polynomial Pr6zbc(x) evaluated at x=0
    // L4(0)'''' = 24
    double d4_c = 24.0 * fx4_c;

    d2dx2 = ((wc/g_c)*(d2_c - g1*d2_1 - g2*d2_2) + w1*d2_1 + w2*d2_2) / (dx*dx);
    
    d4dx4 = ((wc/g_c) * d4_c) / (dx*dx*dx*dx);

}

#endif
