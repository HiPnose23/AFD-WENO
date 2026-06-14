#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>

#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include <AMReX_ParallelDescriptor.H>

#include "../weno_lib.hpp"

// =========================================================
// 1. Setup and Parameters (Made dynamic for the loop)
// =========================================================
double cx = 1.0, cy = 1.0, cz = 1.0;
double lam = 0.0; // Set to 0.0 to test pure advection spatial order
int N;
double dx, dy, dz;

inline double get_x(int i) { return -1.0 + (i + 0.5) * dx; }
inline double get_y(int j) { return -1.0 + (j + 0.5) * dy; }
inline double get_z(int k) { return -1.0 + (k + 0.5) * dz; }

// Global Pointers for AMReX states (Exactly as you had them)
amrex::MultiFab* E_center; 
amrex::MultiFab* E_face;   

void initialize_equilibrium_arrays(const amrex::Geometry& geom) {
    const double pi = std::acos(-1.0);

    for (amrex::MFIter mfi(*E_center); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox(); 
        auto const& Ec  = E_center->array(mfi);
        auto const& Ef  = E_face->array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double x = get_x(i), y = get_y(j), z = get_z(k);
            Ec(i,j,k) = std::exp(lam * (std::sin(pi * x) + std::sin(pi * y) + std::sin(pi * z)));
            
            double x_f = x + 0.5 * dx;
            Ef(i,j,k,0) = std::exp(lam * (std::sin(pi * x_f) + std::sin(pi * y) + std::sin(pi * z)));
            
            double y_f = y + 0.5 * dy;
            Ef(i,j,k,1) = std::exp(lam * (std::sin(pi * x) + std::sin(pi * y_f) + std::sin(pi * z)));
            
            double z_f = z + 0.5 * dz;
            Ef(i,j,k,2) = std::exp(lam * (std::sin(pi * x) + std::sin(pi * y) + std::sin(pi * z_f)));
        });
    }
    E_center->FillBoundary(geom.periodicity());
    E_face->FillBoundary(geom.periodicity());  
}

// =========================================================
// 2. Unified A-WENO Flux Function (Untouched!)
// =========================================================
void compute_wb_flux_3D(const amrex::MultiFab& u, amrex::MultiFab& flux, const amrex::Geometry& geom) {
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox(); 
        
        auto const& u_arr = u.const_array(mfi);
        auto const& Ec  = E_center->const_array(mfi);
        auto const& Ef  = E_face->const_array(mfi); 
        auto const& flx = flux.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double w_L, w_R, dwdx;

            // --- X-Direction ---
            weno_ao_3_interpolation(
                u_arr(i-2,j,k)/Ec(i-2,j,k), u_arr(i-1,j,k)/Ec(i-1,j,k),
                u_arr(i,j,k)/Ec(i,j,k), u_arr(i+1,j,k)/Ec(i+1,j,k), u_arr(i+2,j,k)/Ec(i+2,j,k),
                dx, w_L, w_R, dwdx);
            
            double h_x = cx * w_R * Ef(i,j,k,0);
            double d2f_x = (- (5.0/48.0)*cx*u_arr(i-2,j,k) + (13.0/16.0)*cx*u_arr(i-1,j,k)
                            - (17.0/24.0)*cx*u_arr(i,j,k) - (17.0/24.0)*cx*u_arr(i+1,j,k)
                            + (13.0/16.0)*cx*u_arr(i+2,j,k) - (5.0/48.0)*cx*u_arr(i+3,j,k)) / (dx*dx);
            flx(i,j,k,0) = h_x - (dx*dx/24.0) * d2f_x;

            // --- Y-Direction ---
            weno_ao_3_interpolation(
                u_arr(i,j-2,k)/Ec(i,j-2,k), u_arr(i,j-1,k)/Ec(i,j-1,k),
                u_arr(i,j,k)/Ec(i,j,k), u_arr(i,j+1,k)/Ec(i,j+1,k), u_arr(i,j+2,k)/Ec(i,j+2,k),
                dy, w_L, w_R, dwdx);
            
            double h_y = cy * w_R * Ef(i,j,k,1);
            double d2f_y = (- (5.0/48.0)*cy*u_arr(i,j-2,k) + (13.0/16.0)*cy*u_arr(i,j-1,k)
                            - (17.0/24.0)*cy*u_arr(i,j,k) - (17.0/24.0)*cy*u_arr(i,j+1,k)
                            + (13.0/16.0)*cy*u_arr(i,j+2,k) - (5.0/48.0)*cy*u_arr(i,j+3,k)) / (dy*dy);
            flx(i,j,k,1) = h_y - (dy*dy/24.0) * d2f_y;

            // --- Z-Direction ---
            weno_ao_3_interpolation(
                u_arr(i,j,k-2)/Ec(i,j,k-2), u_arr(i,j,k-1)/Ec(i,j,k-1),
                u_arr(i,j,k)/Ec(i,j,k), u_arr(i,j,k+1)/Ec(i,j,k+1), u_arr(i,j,k+2)/Ec(i,j,k+2),
                dz, w_L, w_R, dwdx);
            
            double h_z = cz * w_R * Ef(i,j,k,2);
            double d2f_z = (- (5.0/48.0)*cz*u_arr(i,j,k-2) + (13.0/16.0)*cz*u_arr(i,j,k-1)
                            - (17.0/24.0)*cz*u_arr(i,j,k) - (17.0/24.0)*cz*u_arr(i,j,k+1)
                            + (13.0/16.0)*cz*u_arr(i,j,k+2) - (5.0/48.0)*cz*u_arr(i,j,k+3)) / (dz*dz);
            flx(i,j,k,2) = h_z - (dz*dz/24.0) * d2f_z;
        });
    }
    flux.FillBoundary(geom.periodicity());
}

// =========================================================
// 3. RHS Assembly (Untouched!)
// =========================================================
void get_rhs_local_equilibrium_3D(const amrex::MultiFab& u, amrex::MultiFab& dudt, const amrex::Geometry& geom) {
    amrex::MultiFab flux_act(u.boxArray(), u.DistributionMap(), 3, 1);
    compute_wb_flux_3D(u, flux_act, geom);

    amrex::MultiFab flux_eq(u.boxArray(), u.DistributionMap(), 3, 1);
    compute_wb_flux_3D(*E_center, flux_eq, geom);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        
        auto const& f_act = flux_act.const_array(mfi); 
        auto const& f_eq  = flux_eq.const_array(mfi);
        auto const& u_arr = u.const_array(mfi);
        auto const& Ec    = E_center->const_array(mfi);
        auto const& rhs   = dudt.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double dFdx_act = (f_act(i,j,k,0) - f_act(i-1,j,k,0)) / dx;
            double dFdy_act = (f_act(i,j,k,1) - f_act(i,j-1,k,1)) / dy;
            double dFdz_act = (f_act(i,j,k,2) - f_act(i,j,k-1,2)) / dz;

            double Div_E_x = (f_eq(i,j,k,0) - f_eq(i-1,j,k,0)) / dx;
            double Div_E_y = (f_eq(i,j,k,1) - f_eq(i,j-1,k,1)) / dy;
            double Div_E_z = (f_eq(i,j,k,2) - f_eq(i,j,k-1,2)) / dz;

            double scaling_fraction = u_arr(i,j,k) / Ec(i,j,k);

            double Sx = scaling_fraction * Div_E_x;
            double Sy = scaling_fraction * Div_E_y;
            double Sz = scaling_fraction * Div_E_z;

            rhs(i,j,k) = -(dFdx_act + dFdy_act + dFdz_act) + (Sx + Sy + Sz);
        });
    }
}

// =========================================================
// 4. Main Program: Looping Convergence Test
// =========================================================
int main(int argc, char* argv[]) {
    amrex::Initialize(argc, argv);
    { 
        amrex::Print() << "Starting Spatial Convergence Test...\n\n";
        
        std::vector<int> N_vals = {16, 32, 64};
        std::vector<double> l1_err(3), linf_err(3);

        for (int step_n = 0; step_n < 3; ++step_n) {
            N = N_vals[step_n];
            dx = 2.0 / N; dy = 2.0 / N; dz = 2.0 / N;
            
            amrex::Print() << "Running grid N = " << N << "...\n";

            amrex::Box domain(amrex::IntVect(0,0,0), amrex::IntVect(N-1,N-1,N-1));
            amrex::RealBox real_box({AMREX_D_DECL(-1.0, -1.0, -1.0)}, {AMREX_D_DECL(1.0, 1.0, 1.0)});
            amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};
            amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

            amrex::BoxArray ba(domain);
            ba.maxSize(32); 
            amrex::DistributionMapping dm(ba);

            E_center = new amrex::MultiFab(ba, dm, 1, 3);
            E_face   = new amrex::MultiFab(ba, dm, 3, 3); 
            initialize_equilibrium_arrays(geom);

            amrex::MultiFab u(ba, dm, 1, 3);
            amrex::MultiFab u1(ba, dm, 1, 3);
            amrex::MultiFab u2(ba, dm, 1, 3);
            amrex::MultiFab rhs(ba, dm, 1, 0);
            amrex::MultiFab u_exact(ba, dm, 1, 0);

            // Initialize a moving sine wave
            const double pi = std::acos(-1.0);
            for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
                const amrex::Box& bx = mfi.tilebox();
                auto const& u_arr = u.array(mfi);
                auto const& u_ex = u_exact.array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    double x = get_x(i), y = get_y(j), z = get_z(k);
                    double wave = std::sin(pi * x) + std::sin(pi * y) + std::sin(pi * z);
                    u_arr(i,j,k) = wave;
                    u_ex(i,j,k)  = wave; // Exact solution after 1 full period is identical
                });
            }
            u.FillBoundary(geom.periodicity()); 

            double t = 0.0;
            double t_end = 2.0; // Advect for 1 full period (length 2.0, speed 1.0)
            double CFL = 0.1;   // Strict CFL so 3rd-order time error doesn't pollute 4th-order space
            double c_speed = std::max({std::abs(cx), std::abs(cy), std::abs(cz)});

            while (t < t_end) {
                double dt = CFL * dx / c_speed;
                if (t + dt > t_end) dt = t_end - t;

                get_rhs_local_equilibrium_3D(u, rhs, geom);
                amrex::MultiFab::Copy(u1, u, 0, 0, 1, 0);                 
                amrex::MultiFab::Saxpy(u1, dt, rhs, 0, 0, 1, 0);          
                u1.FillBoundary(geom.periodicity());

                get_rhs_local_equilibrium_3D(u1, rhs, geom);
                amrex::MultiFab::LinComb(u2, 0.75, u, 0, 0.25, u1, 0, 0, 1, 0); 
                amrex::MultiFab::Saxpy(u2, 0.25 * dt, rhs, 0, 0, 1, 0);
                u2.FillBoundary(geom.periodicity());

                get_rhs_local_equilibrium_3D(u2, rhs, geom);
                amrex::MultiFab::LinComb(u, 1.0/3.0, u, 0, 2.0/3.0, u2, 0, 0, 1, 0);
                amrex::MultiFab::Saxpy(u, (2.0/3.0) * dt, rhs, 0, 0, 1, 0);
                u.FillBoundary(geom.periodicity());

                t += dt;
            }

            // Calculate Error
            amrex::MultiFab error(ba, dm, 1, 0);
            amrex::MultiFab::Copy(error, u, 0, 0, 1, 0);
            amrex::MultiFab::Subtract(error, u_exact, 0, 0, 1, 0);

            double vol = dx * dy * dz;
            linf_err[step_n] = error.norm0(0);
            l1_err[step_n]   = error.norm1(0) * vol;

            delete E_center; delete E_face;
        }

        // Print Convergence Table
        amrex::Print() << "\n=======================================================\n";
        amrex::Print() << "            SPATIAL CONVERGENCE TABLE                  \n";
        amrex::Print() << "=======================================================\n";
        amrex::Print() << std::setw(6) << "N" << std::setw(15) << "L1 Error" << std::setw(12) << "Order" 
                       << std::setw(15) << "L-inf Error" << std::setw(12) << "Order\n";
        amrex::Print() << "-------------------------------------------------------\n";

        for (int i = 0; i < 3; ++i) {
            double order_l1 = (i == 0) ? 0.0 : std::log(l1_err[i-1] / l1_err[i]) / std::log(2.0);
            double order_linf = (i == 0) ? 0.0 : std::log(linf_err[i-1] / linf_err[i]) / std::log(2.0);

            amrex::Print() << std::setw(6) << N_vals[i] 
                           << std::scientific << std::setprecision(4) 
                           << std::setw(15) << l1_err[i]
                           << std::fixed << std::setprecision(2) 
                           << std::setw(12) << (i == 0 ? "-" : std::to_string(order_l1))
                           << std::scientific << std::setprecision(4) 
                           << std::setw(15) << linf_err[i]
                           << std::setw(12) << (i == 0 ? "-" : std::to_string(order_linf)) << "\n";
        }
        amrex::Print() << "=======================================================\n";
    }
    amrex::Finalize();
    return 0;
}
