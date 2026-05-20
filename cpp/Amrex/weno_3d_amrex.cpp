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
// 1. Setup and Parameters
// =========================================================
constexpr double cx = 1.0, cy = 1.0, cz = 1.0;
constexpr double lam = -0.3; 
constexpr int N = 16;

const double dx = 2.0 / N, dy = 2.0 / N, dz = 2.0 / N;

inline double get_x(int i) { return -1.0 + (i + 0.5) * dx; }
inline double get_y(int j) { return -1.0 + (j + 0.5) * dy; }
inline double get_z(int k) { return -1.0 + (k + 0.5) * dz; }

// Global Pointers for AMReX states 
amrex::MultiFab* E_center; // ncomp = 1
amrex::MultiFab* E_face;   // ncomp = 3 (0=X, 1=Y, 2=Z)

void initialize_equilibrium_arrays(const amrex::Geometry& geom) {
    const double pi = std::acos(-1.0);

    for (amrex::MFIter mfi(*E_center); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox(); 
        
        auto const& Ec  = E_center->array(mfi);
        auto const& Ef  = E_face->array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double x = get_x(i), y = get_y(j), z = get_z(k);
            
            // Center Value (Periodic Sine Wave Equilibrium)
            Ec(i,j,k) = std::exp(lam * (std::sin(pi * x) + std::sin(pi * y) + std::sin(pi * z)));
            
            // X-Face (n=0)
            double x_f = x + 0.5 * dx;
            Ef(i,j,k,0) = std::exp(lam * (std::sin(pi * x_f) + std::sin(pi * y) + std::sin(pi * z)));
            
            // Y-Face (n=1)
            double y_f = y + 0.5 * dy;
            Ef(i,j,k,1) = std::exp(lam * (std::sin(pi * x) + std::sin(pi * y_f) + std::sin(pi * z)));
            
            // Z-Face (n=2)
            double z_f = z + 0.5 * dz;
            Ef(i,j,k,2) = std::exp(lam * (std::sin(pi * x) + std::sin(pi * y) + std::sin(pi * z_f)));
        });
    }
    
    // Fill the periodic ghost cells natively via AMReX
    E_center->FillBoundary(geom.periodicity());
    E_face->FillBoundary(geom.periodicity());  
}

// =========================================================
// 2. Unified A-WENO Flux Function (Eq. 2.6 & 2.7 from Paper)
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

            // --- X-Direction Face Flux ---
            weno_ao_3_interpolation(
                u_arr(i-2,j,k)/Ec(i-2,j,k), u_arr(i-1,j,k)/Ec(i-1,j,k),
                u_arr(i,j,k)/Ec(i,j,k), u_arr(i+1,j,k)/Ec(i+1,j,k), u_arr(i+2,j,k)/Ec(i+2,j,k),
                dx, w_L, w_R, dwdx);
            
            double h_x = cx * w_R * Ef(i,j,k,0); // Physical Upwind Flux
            
            // 6-point high-order correction applied directly to the cell-centered physical flux
            double d2f_x = (
                - (5.0/48.0)  * cx * u_arr(i-2,j,k)
                + (13.0/16.0) * cx * u_arr(i-1,j,k)
                - (17.0/24.0) * cx * u_arr(i,j,k)
                - (17.0/24.0) * cx * u_arr(i+1,j,k)
                + (13.0/16.0) * cx * u_arr(i+2,j,k)
                - (5.0/48.0)  * cx * u_arr(i+3,j,k)) / (dx*dx);
            
            flx(i,j,k,0) = h_x - (dx*dx/24.0) * d2f_x;

            // --- Y-Direction Face Flux ---
            weno_ao_3_interpolation(
                u_arr(i,j-2,k)/Ec(i,j-2,k), u_arr(i,j-1,k)/Ec(i,j-1,k),
                u_arr(i,j,k)/Ec(i,j,k), u_arr(i,j+1,k)/Ec(i,j+1,k), u_arr(i,j+2,k)/Ec(i,j+2,k),
                dy, w_L, w_R, dwdx);
            
            double h_y = cy * w_R * Ef(i,j,k,1);
            double d2f_y = (
                - (5.0/48.0)  * cy * u_arr(i,j-2,k)
                + (13.0/16.0) * cy * u_arr(i,j-1,k)
                - (17.0/24.0) * cy * u_arr(i,j,k)
                - (17.0/24.0) * cy * u_arr(i,j+1,k)
                + (13.0/16.0) * cy * u_arr(i,j+2,k)
                - (5.0/48.0)  * cy * u_arr(i,j+3,k)) / (dy*dy);
                
            flx(i,j,k,1) = h_y - (dy*dy/24.0) * d2f_y;

            // --- Z-Direction Face Flux ---
            weno_ao_3_interpolation(
                u_arr(i,j,k-2)/Ec(i,j,k-2), u_arr(i,j,k-1)/Ec(i,j,k-1),
                u_arr(i,j,k)/Ec(i,j,k), u_arr(i,j,k+1)/Ec(i,j,k+1), u_arr(i,j,k+2)/Ec(i,j,k+2),
                dz, w_L, w_R, dwdx);
            
            double h_z = cz * w_R * Ef(i,j,k,2);
            double d2f_z = (
                - (5.0/48.0)  * cz * u_arr(i,j,k-2)
                + (13.0/16.0) * cz * u_arr(i,j,k-1)
                - (17.0/24.0) * cz * u_arr(i,j,k)
                - (17.0/24.0) * cz * u_arr(i,j,k+1)
                + (13.0/16.0) * cz * u_arr(i,j,k+2)
                - (5.0/48.0)  * cz * u_arr(i,j,k+3)) / (dz*dz);
                
            flx(i,j,k,2) = h_z - (dz*dz/24.0) * d2f_z;
        });
    }
    
    // Fill boundaries so divergence calculation can access cell i-1
    flux.FillBoundary(geom.periodicity());
}

// =========================================================
// 3. RHS Assembly with Generalized Source Scaling 
// =========================================================
void get_rhs_local_equilibrium_3D(const amrex::MultiFab& u, amrex::MultiFab& Rhs, const amrex::Geometry& geom) {
    // 1. Calculate Actual Numerical Fluxes: F_num(u)
    amrex::MultiFab flux_act(u.boxArray(), u.DistributionMap(), 3, 1);
    compute_wb_flux_3D(u, flux_act, geom);

    // 2. Calculate Reference Equilibrium Numerical Fluxes: F_num(u_{e,loc})
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
            
            // Actual Flux Divergence
            double dFdx_act = (f_act(i,j,k,0) - f_act(i-1,j,k,0)) / dx;
            double dFdy_act = (f_act(i,j,k,1) - f_act(i,j-1,k,1)) / dy;
            double dFdz_act = (f_act(i,j,k,2) - f_act(i,j,k-1,2)) / dz;

            // Numerical Divergence of the Unscaled Reference State
            double Div_E_x = (f_eq(i,j,k,0) - f_eq(i-1,j,k,0)) / dx;
            double Div_E_y = (f_eq(i,j,k,1) - f_eq(i,j-1,k,1)) / dy;
            double Div_E_z = (f_eq(i,j,k,2) - f_eq(i,j,k-1,2)) / dz;

            // -------------------------------------------------------------
            // Explicit Source Term Scaling Fraction: s(u, x) / s(u_e, x)
            // For linear advection with s(u, x) = Div(cu), the ratio 
            // simplifies algebraically to u_i / E(x_i), ensuring safety 
            // against division by zero at the roots of the sine wave.
            // -------------------------------------------------------------
            double scaling_fraction = u_arr(i,j,k) / Ec(i,j,k);

            // S_i = [ s(u)/s(E) ] * Div( F_num(E) )
            double Sx = scaling_fraction * Div_E_x;
            double Sy = scaling_fraction * Div_E_y;
            double Sz = scaling_fraction * Div_E_z;

            // RHS = -Div(F_act) + S_i
            rhs(i,j,k) = -(dFdx_act + dFdy_act + dFdz_act) + (Sx + Sy + Sz);
        });
    }
}

// =========================================================
// 4. Main Simulation & SSP-RK3 Time Stepping
// =========================================================
int main(int argc, char* argv[]) {
    amrex::Initialize(argc, argv);
    { 
        amrex::Print() << "Initializing AMReX A-WENO 3D SSP-RK3 Solver...\n";
        
        amrex::Box domain(amrex::IntVect(0,0,0), amrex::IntVect(N-1,N-1,N-1));
        amrex::RealBox real_box({AMREX_D_DECL(-1.0, -1.0, -1.0)}, {AMREX_D_DECL(1.0, 1.0, 1.0)});
        amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};
        amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

        amrex::BoxArray ba(domain);
        ba.maxSize(8); 
        amrex::DistributionMapping dm(ba);

        // Allocate equilibrium globals
        E_center = new amrex::MultiFab(ba, dm, 1, 3);
        E_face   = new amrex::MultiFab(ba, dm, 3, 3); 

        initialize_equilibrium_arrays(geom);

        // Allocate State and RK3 arrays
        amrex::MultiFab u(ba, dm, 1, 3);
        amrex::MultiFab u1(ba, dm, 1, 3);
        amrex::MultiFab u2(ba, dm, 1, 3);
        amrex::MultiFab rhs(ba, dm, 1, 0);

        // Set Initial Condition
        amrex::MultiFab::Copy(u, *E_center, 0, 0, 1, 3); 
        u.FillBoundary(geom.periodicity()); 

        double t = 0.0;
        double t_end = 1.0;
        double CFL = 0.4;
        double c_speed = std::max({std::abs(cx), std::abs(cy), std::abs(cz)});

        amrex::Print() << "Starting Time Integration...\n";
        int step = 0;

        while (t < t_end) {
            double dt = CFL * dx / c_speed;
            if (t + dt > t_end) dt = t_end - t;

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
            step++;
            
            if (step % 10 == 0) {
                amrex::Print() << "Step " << step << " | Time: " << t << " / " << t_end << "\n";
            }
        }

        // =========================================================
        // Calculate Machine Precision Well-Balanced Error
        // =========================================================
        amrex::MultiFab error(ba, dm, 1, 0);
        amrex::MultiFab::Copy(error, u, 0, 0, 1, 0);
        amrex::MultiFab::Subtract(error, *E_center, 0, 0, 1, 0);

        double vol = dx * dy * dz;
        double linf_err = error.norm0(0);
        double l1_err   = error.norm1(0) * vol;
        double l2_err   = error.norm2(0) * std::sqrt(vol);

        amrex::Print() << "\n--- Well-Balanced Verification ---\n"
                       << std::scientific << std::setprecision(8)
                       << "L1 Error:    " << l1_err << "\n"
                       << "L2 Error:    " << l2_err << "\n"
                       << "L-inf Error: " << linf_err << "\n"
                       << "----------------------------------\n";
        
        // =========================================================
        // Write Output to AMReX Plotfile
        // =========================================================
        std::string plot_filename = "plt_final";
        amrex::Vector<std::string> var_names = {"u"}; 
        amrex::WriteSingleLevelPlotfile(plot_filename, u, var_names, geom, t, step);
        amrex::Print() << "Successfully wrote plotfile to directory: " << plot_filename << "\n";

        // Cleanup
        delete E_center; delete E_face;
    }
    amrex::Finalize();
    return 0;
}
