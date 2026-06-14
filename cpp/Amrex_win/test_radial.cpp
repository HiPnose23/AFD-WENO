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
#include <AMReX_PlotFileUtil.H>

#include "../weno_lib.hpp"

// =========================================================
// 1. Setup and Parameters (Cylindrical: r, theta, z)
// =========================================================
constexpr double cr = 1.0, ctheta = 1.0, cz = 1.0;
constexpr double lam = -0.3; // Non-zero to ensure equilibrium scaling is active

// Spacing variables are now global but non-const so they can be modified in the convergence loop
double dr, dtheta, dz_step;

// r goes from 1.0 to 3.0 to avoid negative radius
inline double get_r(int i)     { return 1.0 + (i + 0.5) * dr; } 
inline double get_theta(int j) { return -1.0 + (j + 0.5) * dtheta; }
inline double get_z(int k)     { return -1.0 + (k + 0.5) * dz_step; }

// Global Pointers for AMReX states 
amrex::MultiFab* E_center; // ncomp = 1
amrex::MultiFab* E_face;   // ncomp = 3 (0=R, 1=Theta, 2=Z)

// =========================================================
// 2. Analytical Exact Solution (Method of Characteristics)
// =========================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
double get_exact_u(double r, double theta, double z, double t, 
                   double c_r, double c_theta, double c_z, double lambda) {
    
    // Trace back the characteristics to t=0
    double r0 = r - c_r * t;
    double z0 = z - c_z * t;
    double theta0;
    
    if (std::abs(c_r) > 1e-12) {
        theta0 = theta - (c_theta / c_r) * std::log(r / r0);
    } else {
        theta0 = theta - (c_theta / r) * t;
    }
    
    // Initialize a Gaussian pulse centered cleanly inside the domain at (2.0, 0.0, 0.0)
    // using a sharp drop-off (30.0) so it doesn't wrap over the periodic boundaries
    double r_dist = r0 - 2.0;
    double theta_dist = theta0 - 0.0;
    double z_dist = z0 - 0.0;
    double dist2 = r_dist*r_dist + theta_dist*theta_dist + z_dist*z_dist;
    
    double q0 = 1.0 + 0.1 * std::exp(-30.0 * dist2);
    
    double pi = std::acos(-1.0);
    double E_curr = std::exp(lambda * (std::sin(pi * r) + std::sin(pi * theta) + std::sin(pi * z)));
    
    return E_curr * q0;
}

void initialize_equilibrium_arrays(const amrex::Geometry& geom) {
    const double pi = std::acos(-1.0);

    for (amrex::MFIter mfi(*E_center); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox(); 
        
        auto const& Ec  = E_center->array(mfi);
        auto const& Ef  = E_face->array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double r = get_r(i), theta = get_theta(j), z = get_z(k);
            
            // Center Value
            Ec(i,j,k) = std::exp(lam * (std::sin(pi * r) + std::sin(pi * theta) + std::sin(pi * z)));
            
            // R-Face (n=0)
            double r_f = r + 0.5 * dr;
            Ef(i,j,k,0) = std::exp(lam * (std::sin(pi * r_f) + std::sin(pi * theta) + std::sin(pi * z)));
            
            // Theta-Face (n=1)
            double theta_f = theta + 0.5 * dtheta;
            Ef(i,j,k,1) = std::exp(lam * (std::sin(pi * r) + std::sin(pi * theta_f) + std::sin(pi * z)));
            
            // Z-Face (n=2)
            double z_f = z + 0.5 * dz_step;
            Ef(i,j,k,2) = std::exp(lam * (std::sin(pi * r) + std::sin(pi * theta) + std::sin(pi * z_f)));
        });
    }
    
    E_center->FillBoundary(geom.periodicity());
    E_face->FillBoundary(geom.periodicity());  
}

// =========================================================
// 3. Unified A-WENO Flux Function (Cylindrical)
// =========================================================
void compute_wb_flux_3D(const amrex::MultiFab& u, amrex::MultiFab& flux, const amrex::Geometry& geom) {
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox(); 
        
        auto const& u_arr = u.const_array(mfi);
        auto const& Ec  = E_center->const_array(mfi);
        auto const& Ef  = E_face->const_array(mfi); 
        auto const& flx = flux.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double w_L, w_R, dwdr;

            // --- R-Direction Face Flux ---
            weno_ao_3_interpolation(
                u_arr(i-2,j,k)/Ec(i-2,j,k), u_arr(i-1,j,k)/Ec(i-1,j,k),
                u_arr(i,j,k)/Ec(i,j,k), u_arr(i+1,j,k)/Ec(i+1,j,k), u_arr(i+2,j,k)/Ec(i+2,j,k),
                dr, w_L, w_R, dwdr);
            
            double h_r = cr * w_R * Ef(i,j,k,0); 
            double r_face = get_r(i) + 0.5 * dr;
            
            double d2f_r = (
                - (5.0/48.0)  * cr * get_r(i-2) * u_arr(i-2,j,k)
                + (13.0/16.0) * cr * get_r(i-1) * u_arr(i-1,j,k)
                - (17.0/24.0) * cr * get_r(i)   * u_arr(i,j,k)
                - (17.0/24.0) * cr * get_r(i+1) * u_arr(i+1,j,k)
                + (13.0/16.0) * cr * get_r(i+2) * u_arr(i+2,j,k)
                - (5.0/48.0)  * cr * get_r(i+3) * u_arr(i+3,j,k)) / (dr*dr);
            
            flx(i,j,k,0) = (r_face * h_r) - (dr*dr/24.0) * d2f_r;

            // --- Theta-Direction Face Flux ---
            weno_ao_3_interpolation(
                u_arr(i,j-2,k)/Ec(i,j-2,k), u_arr(i,j-1,k)/Ec(i,j-1,k),
                u_arr(i,j,k)/Ec(i,j,k), u_arr(i,j+1,k)/Ec(i,j+1,k), u_arr(i,j+2,k)/Ec(i,j+2,k),
                dtheta, w_L, w_R, dwdr);
            
            double h_theta = ctheta * w_R * Ef(i,j,k,1);
            double d2f_theta = (
                - (5.0/48.0)  * ctheta * u_arr(i,j-2,k)
                + (13.0/16.0) * ctheta * u_arr(i,j-1,k)
                - (17.0/24.0) * ctheta * u_arr(i,j,k)
                - (17.0/24.0) * ctheta * u_arr(i,j+1,k)
                + (13.0/16.0) * ctheta * u_arr(i,j+2,k)
                - (5.0/48.0)  * ctheta * u_arr(i,j+3,k)) / (dtheta*dtheta);
                
            flx(i,j,k,1) = h_theta - (dtheta*dtheta/24.0) * d2f_theta;

            // --- Z-Direction Face Flux ---
            weno_ao_3_interpolation(
                u_arr(i,j,k-2)/Ec(i,j,k-2), u_arr(i,j,k-1)/Ec(i,j,k-1),
                u_arr(i,j,k)/Ec(i,j,k), u_arr(i,j,k+1)/Ec(i,j,k+1), u_arr(i,j,k+2)/Ec(i,j,k+2),
                dz_step, w_L, w_R, dwdr);
            
            double h_z = cz * w_R * Ef(i,j,k,2);
            double d2f_z = (
                - (5.0/48.0)  * cz * u_arr(i,j,k-2)
                + (13.0/16.0) * cz * u_arr(i,j,k-1)
                - (17.0/24.0) * cz * u_arr(i,j,k)
                - (17.0/24.0) * cz * u_arr(i,j,k+1)
                + (13.0/16.0) * cz * u_arr(i,j,k+2)
                - (5.0/48.0)  * cz * u_arr(i,j,k+3)) / (dz_step*dz_step);
                
            flx(i,j,k,2) = h_z - (dz_step*dz_step/24.0) * d2f_z;
        });
    }
    
    flux.FillBoundary(geom.periodicity());
}

// =========================================================
// 4. RHS Assembly 
// =========================================================
void get_rhs_local_equilibrium_3D(const amrex::MultiFab& u, amrex::MultiFab& Rhs, const amrex::Geometry& geom) {
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
        auto const& rhs   = Rhs.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double r_i = get_r(i); 

            double dFdr_act = (f_act(i,j,k,0) - f_act(i-1,j,k,0)) / (r_i * dr);     
            double dFdt_act = (f_act(i,j,k,1) - f_act(i,j-1,k,1)) / (r_i * dtheta); 
            double dFdz_act = (f_act(i,j,k,2) - f_act(i,j,k-1,2)) / dz_step;        

            double Div_E_r = (f_eq(i,j,k,0) - f_eq(i-1,j,k,0)) / (r_i * dr);
            double Div_E_t = (f_eq(i,j,k,1) - f_eq(i,j-1,k,1)) / (r_i * dtheta);
            double Div_E_z = (f_eq(i,j,k,2) - f_eq(i,j,k-1,2)) / dz_step;

            double scaling_fraction = u_arr(i,j,k) / Ec(i,j,k);

            double Sr = scaling_fraction * Div_E_r;
            double St = scaling_fraction * Div_E_t;
            double Sz = scaling_fraction * Div_E_z;

            rhs(i,j,k) = -(dFdr_act + dFdt_act + dFdz_act) + (Sr + St + Sz);
        });
    }
}

// =========================================================
// 5. Main Simulation: Grid Convergence Study
// =========================================================
int main(int argc, char* argv[]) {
    amrex::Initialize(argc, argv);
    { 
        amrex::Print() << "Initializing AMReX A-WENO 3D Convergence Test...\n";
        
        std::vector<int> N_vals = {32, 64, 128};
        std::vector<double> l1_errs;
        
        for (int N : N_vals) {
            amrex::Print() << "\n--- Running for N = " << N << " ---\n";
            
            // Reassign grid spacings dynamically
            dr = 2.0 / N; dtheta = 2.0 / N; dz_step = 2.0 / N;
            
            amrex::Box domain(amrex::IntVect(0,0,0), amrex::IntVect(N-1,N-1,N-1));
            amrex::RealBox real_box({AMREX_D_DECL(1.0, -1.0, -1.0)}, {AMREX_D_DECL(3.0, 1.0, 1.0)});
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

            // Set Initial Condition with the Exact Gaussian Perturbation
            for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
                const amrex::Box& bx = mfi.validbox();
                auto const& u_arr = u.array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    double r = get_r(i), theta = get_theta(j), z = get_z(k);
                    u_arr(i,j,k) = get_exact_u(r, theta, z, 0.0, cr, ctheta, cz, lam);
                });
            }
            u.FillBoundary(geom.periodicity()); 

            double t = 0.0;
            double t_end = 0.1; // Short time prevents the pulse from hitting boundaries
            double CFL = 0.4;
            double c_speed = std::max({std::abs(cr), std::abs(ctheta), std::abs(cz)});

            while (t < t_end) {
                double dt = CFL * dr / c_speed;
                // Clamp last dt perfectly to avoid overshooting exact analytical state
                if (t + dt > t_end - 1e-12) dt = t_end - t;

                // --- RK3 Stage 1 ---
                get_rhs_local_equilibrium_3D(u, rhs, geom);
                amrex::MultiFab::Copy(u1, u, 0, 0, 1, 0);                 
                amrex::MultiFab::Saxpy(u1, dt, rhs, 0, 0, 1, 0);          
                u1.FillBoundary(geom.periodicity());

                // --- RK3 Stage 2 ---
                get_rhs_local_equilibrium_3D(u1, rhs, geom);
                amrex::MultiFab::LinComb(u2, 0.75, u, 0, 0.25, u1, 0, 0, 1, 0); 
                amrex::MultiFab::Saxpy(u2, 0.25 * dt, rhs, 0, 0, 1, 0);
                u2.FillBoundary(geom.periodicity());

                // --- RK3 Stage 3 ---
                get_rhs_local_equilibrium_3D(u2, rhs, geom);
                amrex::MultiFab::LinComb(u, 1.0/3.0, u, 0, 2.0/3.0, u2, 0, 0, 1, 0);
                amrex::MultiFab::Saxpy(u, (2.0/3.0) * dt, rhs, 0, 0, 1, 0);
                u.FillBoundary(geom.periodicity());

                t += dt;
            }

            // =========================================================
            // Evaluate Error Against Exact Analytical State
            // =========================================================
            amrex::MultiFab error(ba, dm, 1, 0);
            for (amrex::MFIter mfi(error); mfi.isValid(); ++mfi) {
                const amrex::Box& bx = mfi.validbox();
                auto const& err_arr = error.array(mfi);
                auto const& u_arr = u.const_array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    double r = get_r(i), theta = get_theta(j), z = get_z(k);
                    double u_exact = get_exact_u(r, theta, z, t, cr, ctheta, cz, lam);
                    err_arr(i,j,k) = std::abs(u_arr(i,j,k) - u_exact);
                });
            }

            double vol = dr * dtheta * dz_step; 
            double l1_err = error.norm1(0) * vol;
            l1_errs.push_back(l1_err);

            amrex::Print() << "L1 Error: " << l1_err << "\n";
            
            delete E_center; delete E_face;
        }

        amrex::Print() << "\n=== Order of Accuracy Analysis ===\n";
        for (size_t i = 1; i < N_vals.size(); ++i) {
            double order = std::log2(l1_errs[i-1] / l1_errs[i]);
            amrex::Print() << "Order N=" << N_vals[i-1] << " to N=" << N_vals[i] << " : " << order << "\n";
        }
    }
    amrex::Finalize();
    return 0;
}
