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
constexpr double lam = -0.3; 
constexpr int N = 1024;


const double dr = 2.0 / N, dtheta = 2.0 / N, dz_step = 2.0 / N;

// r goes from 1.0 to 3.0 to avoid negative radius
inline double get_r(int i)     { return 1.0 + (i + 0.5) * dr; } 
inline double get_theta(int j) { return -1.0 + (j + 0.5) * dtheta; }
inline double get_z(int k)     { return -1.0 + (k + 0.5) * dz_step; }

// Global Pointers for AMReX states 
amrex::MultiFab* E_center; // ncomp = 1
amrex::MultiFab* E_face;   // ncomp = 3 (0=R, 1=Theta, 2=Z)

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
    
    // Fill the periodic ghost cells natively via AMReX
    E_center->FillBoundary(geom.periodicity());
    E_face->FillBoundary(geom.periodicity());  
}

// =========================================================
// 2. Unified A-WENO Flux Function (Cylindrical)
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

            // --- NEW: R-Direction Face Flux (Cylindrical Geometry) ---
            weno_ao_3_interpolation(
                u_arr(i-2,j,k)/Ec(i-2,j,k), u_arr(i-1,j,k)/Ec(i-1,j,k),
                u_arr(i,j,k)/Ec(i,j,k), u_arr(i+1,j,k)/Ec(i+1,j,k), u_arr(i+2,j,k)/Ec(i+2,j,k),
                dr, w_L, w_R, dwdr);
            
            double h_r = cr * w_R * Ef(i,j,k,0); // Physical Upwind Flux
            
            // 1. Get the exact radius at the face
            double r_face = get_r(i) + 0.5 * dr;
            
            // 2. 6-point high-order correction applied to the GEOMETRIC flux (r * f)
            // Note: Each point in the stencil is multiplied by its own local radius!
            double d2f_r = (
                - (5.0/48.0)  * cr * get_r(i-2) * u_arr(i-2,j,k)
                + (13.0/16.0) * cr * get_r(i-1) * u_arr(i-1,j,k)
                - (17.0/24.0) * cr * get_r(i)   * u_arr(i,j,k)
                - (17.0/24.0) * cr * get_r(i+1) * u_arr(i+1,j,k)
                + (13.0/16.0) * cr * get_r(i+2) * u_arr(i+2,j,k)
                - (5.0/48.0)  * cr * get_r(i+3) * u_arr(i+3,j,k)) / (dr*dr);
            
            // 3. Final face flux is scaled by the face radius
            flx(i,j,k,0) = (r_face * h_r) - (dr*dr/24.0) * d2f_r;

            // --- Theta-Direction Face Flux (No r metric inside derivative) ---
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

            // --- Z-Direction Face Flux (No r metric inside derivative) ---
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
            
            double r_i = get_r(i); // Cell center radius

            // Actual Flux Divergence
            double dFdr_act = (f_act(i,j,k,0) - f_act(i-1,j,k,0)) / (r_i * dr);     // Radial: (1/r) * d(rF)/dr
            double dFdt_act = (f_act(i,j,k,1) - f_act(i,j-1,k,1)) / (r_i * dtheta); // Azimuthal: (1/r) * dF/dtheta
            double dFdz_act = (f_act(i,j,k,2) - f_act(i,j,k-1,2)) / dz_step;        // Axial: dF/dz (Unchanged)

            // Numerical Divergence of the Unscaled Reference State
            double Div_E_r = (f_eq(i,j,k,0) - f_eq(i-1,j,k,0)) / (r_i * dr);
            double Div_E_t = (f_eq(i,j,k,1) - f_eq(i,j-1,k,1)) / (r_i * dtheta);
            double Div_E_z = (f_eq(i,j,k,2) - f_eq(i,j,k-1,2)) / dz_step;

            // -------------------------------------------------------------
            // Explicit Source Term Scaling Fraction
            // -------------------------------------------------------------
            double scaling_fraction = u_arr(i,j,k) / Ec(i,j,k);

            // S_i = [ s(u)/s(E) ] * Div( F_num(E) )
            double Sr = scaling_fraction * Div_E_r;
            double St = scaling_fraction * Div_E_t;
            double Sz = scaling_fraction * Div_E_z;

            // RHS = -Div(F_act) + S_i
            rhs(i,j,k) = -(dFdr_act + dFdt_act + dFdz_act) + (Sr + St + Sz);
        });
    }
}

// =========================================================
// 4. Main Simulation & SSP-RK3 Time Stepping
// =========================================================
int main(int argc, char* argv[]) {
    amrex::Initialize(argc, argv);
    { 
        amrex::Print() << "Initializing AMReX A-WENO 3D Cylindrical SSP-RK3 Solver...\n";
        
        amrex::Box domain(amrex::IntVect(0,0,0), amrex::IntVect(N-1,N-1,N-1));
        
        // --- NEW: Domain strictly positive for r ---
        // r = [1.0, 3.0], theta = [-1.0, 1.0], z = [-1.0, 1.0]
        amrex::RealBox real_box({AMREX_D_DECL(1.0, -1.0, -1.0)}, {AMREX_D_DECL(3.0, 1.0, 1.0)});
        
        amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};
        amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

        amrex::BoxArray ba(domain);
        ba.maxSize(128); 
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
        double c_speed = std::max({std::abs(cr), std::abs(ctheta), std::abs(cz)});

        amrex::Print() << "Starting Time Integration...\n";
        BL_PROFILE_VAR("Evolution_Loop", pmain);
        int step = 0;

        while (t < t_end) {
            // Determine max safe dt. Notice we account for cell width dx=dr
            double dt = CFL * dr / c_speed;
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
        BL_PROFILE_VAR_STOP(pmain); 


        // =========================================================
        // Calculate Machine Precision Well-Balanced Error
        // =========================================================
        amrex::MultiFab error(ba, dm, 1, 0);
        amrex::MultiFab::Copy(error, u, 0, 0, 1, 0);
        amrex::MultiFab::Subtract(error, *E_center, 0, 0, 1, 0);

        // Vol in cylindrical coordinates is r * dr * dtheta * dz. Let's just use grid volume for the L1/L2 calculation metric
        double vol = dr * dtheta * dz_step; 
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
