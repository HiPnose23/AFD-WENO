#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>

#include "weno_lib.hpp"

// =========================================================
// 1. Setup and Parameters
// =========================================================
constexpr double cx = 1.0, cy = 1.0, cz = 1.0;
constexpr double lam = -0.3; 
constexpr int N = 16;
constexpr int total_cells = N * N * N;

const double dx = 2.0 / N, dy = 2.0 / N, dz = 2.0 / N;

// 3D to 1D mapping with periodic boundary wrap-around
inline int idx(int i, int j, int k) {
    i = (i % N + N) % N; j = (j % N + N) % N; k = (k % N + N) % N;
    return i + N * (j + N * k);
}

inline double get_x(int i) { return -1.0 + (i + 0.5) * dx; }
inline double get_y(int j) { return -1.0 + (j + 0.5) * dy; }
inline double get_z(int k) { return -1.0 + (k + 0.5) * dz; }

// Global Arrays for equilibrium state 
std::vector<double> E_center(total_cells);
std::vector<double> E_face_x(total_cells), E_face_y(total_cells), E_face_z(total_cells);
std::vector<double> E_prime_x(total_cells), E_prime_y(total_cells), E_prime_z(total_cells);

void initialize_equilibrium_arrays() {
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                int c = idx(i,j,k);
                double x = get_x(i), y = get_y(j), z = get_z(k);
                
                E_center[c] = std::exp((lam / 3.0) * (x/cx + y/cy + z/cz));
                
                E_face_x[c] = std::exp((lam / 3.0) * ((x+0.5*dx)/cx + y/cy + z/cz));
                E_face_y[c] = std::exp((lam / 3.0) * (x/cx + (y+0.5*dy)/cy + z/cz));
                E_face_z[c] = std::exp((lam / 3.0) * (x/cx + y/cy + (z+0.5*dz)/cz));
                
                E_prime_x[c] = (lam / (3.0 * cx)) * E_center[c];
                E_prime_y[c] = (lam / (3.0 * cy)) * E_center[c];
                E_prime_z[c] = (lam / (3.0 * cz)) * E_center[c];
            }
        }
    }
}

// =========================================================
// 2. Well-Balanced 3D Flux Function
// =========================================================
void compute_wb_flux_3D(const std::vector<double>& u, 
                        std::vector<double>& fx, 
                        std::vector<double>& fy, 
                        std::vector<double>& fz) 
{
    std::vector<double> F_star_x(total_cells), fc_x(total_cells);
    std::vector<double> F_star_y(total_cells), fc_y(total_cells);
    std::vector<double> F_star_z(total_cells), fc_z(total_cells);

    // WENO-AO(3) Interpolation
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                int c = idx(i,j,k);
                double w_L, w_R, dwdx;

                // X-Direction
                weno_ao_3_interpolation(
                    u[idx(i-2,j,k)]/E_center[idx(i-2,j,k)], u[idx(i-1,j,k)]/E_center[idx(i-1,j,k)],
                    u[c]/E_center[c], u[idx(i+1,j,k)]/E_center[idx(i+1,j,k)], u[idx(i+2,j,k)]/E_center[idx(i+2,j,k)],
                    dx, w_L, w_R, dwdx);
                F_star_x[c] = cx * w_R * E_face_x[c];
                fc_x[c]     = cx * (dwdx * E_center[c] + (u[c]/E_center[c]) * E_prime_x[c]);

                // Y-Direction
                weno_ao_3_interpolation(
                    u[idx(i,j-2,k)]/E_center[idx(i,j-2,k)], u[idx(i,j-1,k)]/E_center[idx(i,j-1,k)],
                    u[c]/E_center[c], u[idx(i,j+1,k)]/E_center[idx(i,j+1,k)], u[idx(i,j+2,k)]/E_center[idx(i,j+2,k)],
                    dy, w_L, w_R, dwdx);
                F_star_y[c] = cy * w_R * E_face_y[c];
                fc_y[c]     = cy * (dwdx * E_center[c] + (u[c]/E_center[c]) * E_prime_y[c]);

                // Z-Direction
                weno_ao_3_interpolation(
                    u[idx(i,j,k-2)]/E_center[idx(i,j,k-2)], u[idx(i,j,k-1)]/E_center[idx(i,j,k-1)],
                    u[c]/E_center[c], u[idx(i,j,k+1)]/E_center[idx(i,j,k+1)], u[idx(i,j,k+2)]/E_center[idx(i,j,k+2)],
                    dz, w_L, w_R, dwdx);
                F_star_z[c] = cz * w_R * E_face_z[c];
                fc_z[c]     = cz * (dwdx * E_center[c] + (u[c]/E_center[c]) * E_prime_z[c]);
            }
        }
    }

    // Boundary Correction
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                int c = idx(i,j,k);
                double d1dx, d3dx;

                weno_ao_43_boundary(fc_x[idx(i-1,j,k)], fc_x[c], fc_x[idx(i+1,j,k)], fc_x[idx(i+2,j,k)], dx, d1dx, d3dx);
                fx[c] = F_star_x[c] - (dx*dx/24.0) * d1dx;

                weno_ao_43_boundary(fc_y[idx(i,j-1,k)], fc_y[c], fc_y[idx(i,j+1,k)], fc_y[idx(i,j+2,k)], dy, d1dx, d3dx);
                fy[c] = F_star_y[c] - (dy*dy/24.0) * d1dx;

                weno_ao_43_boundary(fc_z[idx(i,j,k-1)], fc_z[c], fc_z[idx(i,j,k+1)], fc_z[idx(i,j,k+2)], dz, d1dx, d3dx);
                fz[c] = F_star_z[c] - (dz*dz/24.0) * d1dx;
            }
        }
    }
}

// =========================================================
// 3. Local Equilibrium RHS 
// =========================================================
std::vector<double> get_rhs_local_equilibrium_3D(const std::vector<double>& u) 
{
    std::vector<double> dudt(total_cells, 0.0);
    
    //  Actual Physical Fluxes
    std::vector<double> fx_act(total_cells), fy_act(total_cells), fz_act(total_cells);
    compute_wb_flux_3D(u, fx_act, fy_act, fz_act);

    // Arrays allocated once for the imaginary loop
    std::vector<double> u_eq_local(total_cells);
    std::vector<double> fx_im(total_cells), fy_im(total_cells), fz_im(total_cells);

    //  Cell-by-cell Local Equilibrium Evaluation
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                int c = idx(i,j,k);
                
                // Actual divergence
                double dFdx_act = (fx_act[c] - fx_act[idx(i-1,j,k)]) / dx;
                double dFdy_act = (fy_act[c] - fy_act[idx(i,j-1,k)]) / dy;
                double dFdz_act = (fz_act[c] - fz_act[idx(i,j,k-1)]) / dz;

                // ----------------------------------------------------
                // TRUE A-WENO SOURCE TERM (Skipping WENO!)
                // ----------------------------------------------------
                double w_j = u[c] / E_center[c];
                
                // We evaluate the algebraic stencils directly 
                // at the right face (i+1/2) and left face (i-1/2)
                // X direction 
                // Right face (i+1/2)
                double f_eq_x_R = cx * w_j * E_face_x[c]; // Trivial Riemann solver
                double d2f_x_R = (
                    - (5.0/48.0)  * cx * w_j * E_center[idx(i-2,j,k)]
                    + (13.0/16.0) * cx * w_j * E_center[idx(i-1,j,k)]
                    - (17.0/24.0) * cx * w_j * E_center[c]
                    - (17.0/24.0) * cx * w_j * E_center[idx(i+1,j,k)]
                    + (13.0/16.0) * cx * w_j * E_center[idx(i+2,j,k)]
                    - (5.0/48.0)  * cx * w_j * E_center[idx(i+3,j,k)]
                ) / (dx*dx);
                double F_eq_x_R = f_eq_x_R - (dx*dx / 24.0) * d2f_x_R;

                // Left face (i-1/2)
                double f_eq_x_L = cx * w_j * E_face_x[idx(i-1,j,k)]; 
                double d2f_x_L = (
                    - (5.0/48.0)  * cx * w_j * E_center[idx(i-3,j,k)]
                    + (13.0/16.0) * cx * w_j * E_center[idx(i-2,j,k)]
                    - (17.0/24.0) * cx * w_j * E_center[idx(i-1,j,k)]
                    - (17.0/24.0) * cx * w_j * E_center[c]
                    + (13.0/16.0) * cx * w_j * E_center[idx(i+1,j,k)]
                    - (5.0/48.0)  * cx * w_j * E_center[idx(i+2,j,k)]
                ) / (dx*dx);
                double F_eq_x_L = f_eq_x_L - (dx*dx / 24.0) * d2f_x_L;
                
                double Sx = (F_eq_x_R - F_eq_x_L) / dx;

                // (Repeat identical algebraic logic for Sy and Sz)
                // Y direction
                // Right face (j+1/2)
                double f_eq_y_R = cy * w_j * E_face_y[c]; 
                double d2f_y_R = (
                    - (5.0/48.0)  * cy * w_j * E_center[idx(i,j-2,k)]
                    + (13.0/16.0) * cy * w_j * E_center[idx(i,j-1,k)]
                    - (17.0/24.0) * cy * w_j * E_center[c]
                    - (17.0/24.0) * cy * w_j * E_center[idx(i,j+1,k)]
                    + (13.0/16.0) * cy * w_j * E_center[idx(i,j+2,k)]
                    - (5.0/48.0)  * cy * w_j * E_center[idx(i,j+3,k)]
                ) / (dy*dy);
                double F_eq_y_R = f_eq_y_R - (dy*dy / 24.0) * d2f_y_R;

                // Left face (j-1/2)
                double f_eq_y_L = cy * w_j * E_face_y[idx(i,j-1,k)]; 
                double d2f_y_L = (
                    - (5.0/48.0)  * cy * w_j * E_center[idx(i,j-3,k)]
                    + (13.0/16.0) * cy * w_j * E_center[idx(i,j-2,k)]
                    - (17.0/24.0) * cy * w_j * E_center[idx(i,j-1,k)]
                    - (17.0/24.0) * cy * w_j * E_center[c]
                    + (13.0/16.0) * cy * w_j * E_center[idx(i,j+1,k)]
                    - (5.0/48.0)  * cy * w_j * E_center[idx(i,j+2,k)]
                ) / (dy*dy);
                double F_eq_y_L = f_eq_y_L - (dy*dy / 24.0) * d2f_y_L;
                
                double Sy = (F_eq_y_R - F_eq_y_L) / dy;

                // Z direction
                // Right face (k+1/2)
                double f_eq_z_R = cz * w_j * E_face_z[c]; 
                double d2f_z_R = (
                    - (5.0/48.0)  * cz * w_j * E_center[idx(i,j,k-2)]
                    + (13.0/16.0) * cz * w_j * E_center[idx(i,j,k-1)]
                    - (17.0/24.0) * cz * w_j * E_center[c]
                    - (17.0/24.0) * cz * w_j * E_center[idx(i,j,k+1)]
                    + (13.0/16.0) * cz * w_j * E_center[idx(i,j,k+2)]
                    - (5.0/48.0)  * cz * w_j * E_center[idx(i,j,k+3)]
                ) / (dz*dz);
                double F_eq_z_R = f_eq_z_R - (dz*dz / 24.0) * d2f_z_R;

                // Left face (k-1/2)
                double f_eq_z_L = cz * w_j * E_face_z[idx(i,j,k-1)]; 
                double d2f_z_L = (
                    - (5.0/48.0)  * cz * w_j * E_center[idx(i,j,k-3)]
                    + (13.0/16.0) * cz * w_j * E_center[idx(i,j,k-2)]
                    - (17.0/24.0) * cz * w_j * E_center[idx(i,j,k-1)]
                    - (17.0/24.0) * cz * w_j * E_center[c]
                    + (13.0/16.0) * cz * w_j * E_center[idx(i,j,k+1)]
                    - (5.0/48.0)  * cz * w_j * E_center[idx(i,j,k+2)]
                ) / (dz*dz);
                double F_eq_z_L = f_eq_z_L - (dz*dz / 24.0) * d2f_z_L;
                
                double Sz = (F_eq_z_R - F_eq_z_L) / dz;

                dudt[c] = -(dFdx_act + dFdy_act + dFdz_act) + (Sx + Sy + Sz);
            }
        }
    }

    return dudt;
}

// =========================================================
// 4. Verification Test
// =========================================================
int main() 
{
    std::cout << "Initializing Pure C++ 3D Literal Local Equilibrium Setup...\n";
    initialize_equilibrium_arrays();

    // Set initial condition to steady state
    std::vector<double> u = E_center;

    // Evaluate RHS
    std::vector<double> rhs = get_rhs_local_equilibrium_3D(u);

    // Track max error (ignoring periodic boundaries to avoid wrap-around shock)
    double max_err = 0.0;
    for (int k = 3; k < N-3; ++k) {
        for (int j = 3; j < N-3; ++j) {
            for (int i = 3; i < N-3; ++i) {
                max_err = std::max(max_err, std::abs(rhs[idx(i,j,k)]));
            }
        }
    }

    std::cout << std::scientific << std::setprecision(10);
    std::cout << "\n--- 3D EXACT STEADY-STATE TEST ---\n";
    std::cout << "Approach: Local Equilibrium Evaluation (General Xing-Shu)\n";
    std::cout << "Max 3D RHS Error: " << max_err << "\n";
    
    if (max_err < 1e-12) {
        std::cout << "SUCCESS: The scheme is well-balanced to machine precision!\n";
    }

    return 0;
}
