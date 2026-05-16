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
amrex::MultiFab* E_prime;  // ncomp = 3 (0=X, 1=Y, 2=Z)

void initialize_equilibrium_arrays(const amrex::Geometry& geom) {
    for (amrex::MFIter mfi(*E_center); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox(); 
        
        auto const& Ec  = E_center->array(mfi);
        auto const& Ef  = E_face->array(mfi);
        auto const& Ep  = E_prime->array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double x = get_x(i), y = get_y(j), z = get_z(k);
            
            Ec(i,j,k) = std::exp((lam / 3.0) * (x/cx + y/cy + z/cz));
            
            // X-component (n=0)
            Ef(i,j,k,0) = std::exp((lam / 3.0) * ((x+0.5*dx)/cx + y/cy + z/cz));
            Ep(i,j,k,0) = (lam / (3.0 * cx)) * Ec(i,j,k);
            
            // Y-component (n=1)
            Ef(i,j,k,1) = std::exp((lam / 3.0) * (x/cx + (y+0.5*dy)/cy + z/cz));
            Ep(i,j,k,1) = (lam / (3.0 * cy)) * Ec(i,j,k);
            
            // Z-component (n=2)
            Ef(i,j,k,2) = std::exp((lam / 3.0) * (x/cx + y/cy + (z+0.5*dz)/cz));
            Ep(i,j,k,2) = (lam / (3.0 * cz)) * Ec(i,j,k);
        });
    }
    

    E_center->FillBoundary(geom.periodicity());
    E_face->FillBoundary(geom.periodicity());  
    E_prime->FillBoundary(geom.periodicity()); 
}

// =========================================================
// 2. Well-Balanced 3D Flux Function
// =========================================================
void compute_wb_flux_3D(const amrex::MultiFab& u, amrex::MultiFab& flux, const amrex::Geometry& geom) {
    // Intermediate temporary arrays (ncomp = 3)
    amrex::MultiFab F_star(u.boxArray(), u.DistributionMap(), 3, 0); 
    amrex::MultiFab fc(u.boxArray(), u.DistributionMap(), 3, 2);

    //  WENO-AO(3) Interpolation
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        
        auto const& u_arr = u.const_array(mfi);
        auto const& Ec  = E_center->const_array(mfi);
        auto const& Ef  = E_face->const_array(mfi); 
        auto const& Ep  = E_prime->const_array(mfi); 

        auto const& Fs  = F_star.array(mfi); 
        auto const& f_c = fc.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double w_L, w_R, dwdx;

            // X-Direction (n=0)
            weno_ao_3_interpolation(
                u_arr(i-2,j,k)/Ec(i-2,j,k), u_arr(i-1,j,k)/Ec(i-1,j,k),
                u_arr(i,j,k)/Ec(i,j,k), u_arr(i+1,j,k)/Ec(i+1,j,k), u_arr(i+2,j,k)/Ec(i+2,j,k),
                dx, w_L, w_R, dwdx);
            Fs(i,j,k,0)  = cx * w_R * Ef(i,j,k,0);
            f_c(i,j,k,0) = cx * (dwdx * Ec(i,j,k) + (u_arr(i,j,k)/Ec(i,j,k)) * Ep(i,j,k,0));

            // Y-Direction (n=1)
            weno_ao_3_interpolation(
                u_arr(i,j-2,k)/Ec(i,j-2,k), u_arr(i,j-1,k)/Ec(i,j-1,k),
                u_arr(i,j,k)/Ec(i,j,k), u_arr(i,j+1,k)/Ec(i,j+1,k), u_arr(i,j+2,k)/Ec(i,j+2,k),
                dy, w_L, w_R, dwdx);
            Fs(i,j,k,1)  = cy * w_R * Ef(i,j,k,1);
            f_c(i,j,k,1) = cy * (dwdx * Ec(i,j,k) + (u_arr(i,j,k)/Ec(i,j,k)) * Ep(i,j,k,1));

            // Z-Direction (n=2)
            weno_ao_3_interpolation(
                u_arr(i,j,k-2)/Ec(i,j,k-2), u_arr(i,j,k-1)/Ec(i,j,k-1),
                u_arr(i,j,k)/Ec(i,j,k), u_arr(i,j,k+1)/Ec(i,j,k+1), u_arr(i,j,k+2)/Ec(i,j,k+2),
                dz, w_L, w_R, dwdx);
            Fs(i,j,k,2)  = cz * w_R * Ef(i,j,k,2);
            f_c(i,j,k,2) = cz * (dwdx * Ec(i,j,k) + (u_arr(i,j,k)/Ec(i,j,k)) * Ep(i,j,k,2));
        });
    }

    // Pass boundary info for all 3 axes at once
    fc.FillBoundary(geom.periodicity());

    // Boundary Correction
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        
        auto const& f_c = fc.const_array(mfi);
        auto const& Fs  = F_star.const_array(mfi);
        auto const& flx = flux.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double d1dx, d3dx;

            // X-correction
            weno_ao_43_boundary(f_c(i-1,j,k,0), f_c(i,j,k,0), f_c(i+1,j,k,0), f_c(i+2,j,k,0), dx, d1dx, d3dx);
            flx(i,j,k,0) = Fs(i,j,k,0) - (dx*dx/24.0) * d1dx;

            // Y-correction
            weno_ao_43_boundary(f_c(i,j-1,k,1), f_c(i,j,k,1), f_c(i,j+1,k,1), f_c(i,j+2,k,1), dy, d1dx, d3dx);
            flx(i,j,k,1) = Fs(i,j,k,1) - (dy*dy/24.0) * d1dx;

            // Z-correction
            weno_ao_43_boundary(f_c(i,j,k-1,2), f_c(i,j,k,2), f_c(i,j,k+1,2), f_c(i,j,k+2,2), dz, d1dx, d3dx);
            flx(i,j,k,2) = Fs(i,j,k,2) - (dz*dz/24.0) * d1dx;
        });
    }

    // Exchange all physical fluxes simultaneously 
    flux.FillBoundary(geom.periodicity());
}

// =========================================================
// 3. Local Equilibrium RHS 
// =========================================================
void get_rhs_local_equilibrium_3D(const amrex::MultiFab& u, amrex::MultiFab& dudt, const amrex::Geometry& geom) {
    // Actual Physical Fluxes (Allocated with 1 ghost cell for divergence)
    amrex::MultiFab flux_act(u.boxArray(), u.DistributionMap(), 3, 1);
    compute_wb_flux_3D(u, flux_act, geom);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        
        auto const& flx = flux_act.const_array(mfi); 
        auto const& u_arr = u.const_array(mfi);
        auto const& Ec = E_center->const_array(mfi);
        auto const& Ef = E_face->const_array(mfi); 
        
        auto const& rhs = dudt.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            
            // Actual divergence
            double dFdx_act = (flx(i,j,k,0) - flx(i-1,j,k,0)) / dx;
            double dFdy_act = (flx(i,j,k,1) - flx(i,j-1,k,1)) / dy;
            double dFdz_act = (flx(i,j,k,2) - flx(i,j,k-1,2)) / dz;

            double w_j = u_arr(i,j,k) / Ec(i,j,k);
            
            // X direction
            double f_eq_x_R = cx * w_j * Ef(i,j,k,0);
            double d2f_x_R = (
                - (5.0/48.0)  * cx * w_j * Ec(i-2,j,k) 
                + (13.0/16.0) * cx * w_j * Ec(i-1,j,k)
                - (17.0/24.0) * cx * w_j * Ec(i,j,k)
                - (17.0/24.0) * cx * w_j * Ec(i+1,j,k) 
                + (13.0/16.0) * cx * w_j * Ec(i+2,j,k)
                - (5.0/48.0)  * cx * w_j * Ec(i+3,j,k)) 
                / (dx*dx);
            double F_eq_x_R = f_eq_x_R - (dx*dx / 24.0) * d2f_x_R;

            double f_eq_x_L = cx * w_j * Ef(i-1,j,k,0);
            double d2f_x_L = (
                - (5.0/48.0)  * cx * w_j * Ec(i-3,j,k) 
                + (13.0/16.0) * cx * w_j * Ec(i-2,j,k) 
                - (17.0/24.0) * cx * w_j * Ec(i-1,j,k)
                - (17.0/24.0) * cx * w_j * Ec(i,j,k)
                + (13.0/16.0) * cx * w_j * Ec(i+1,j,k) 
                - (5.0/48.0)  * cx * w_j * Ec(i+2,j,k))
                / (dx*dx);
            double F_eq_x_L = f_eq_x_L - (dx*dx / 24.0) * d2f_x_L;
            double Sx = (F_eq_x_R - F_eq_x_L) / dx;

            // Y direction
            double f_eq_y_R = cy * w_j * Ef(i,j,k,1);
            double d2f_y_R = (
                - (5.0/48.0)  * cy * w_j * Ec(i,j-2,k) 
                + (13.0/16.0) * cy * w_j * Ec(i,j-1,k) 
                - (17.0/24.0) * cy * w_j * Ec(i,j,k)
                - (17.0/24.0) * cy * w_j * Ec(i,j+1,k) 
                + (13.0/16.0) * cy * w_j * Ec(i,j+2,k) 
                - (5.0/48.0)  * cy * w_j * Ec(i,j+3,k))
                / (dy*dy);
            double F_eq_y_R = f_eq_y_R - (dy*dy / 24.0) * d2f_y_R;

            double f_eq_y_L = cy * w_j * Ef(i,j-1,k,1);
            double d2f_y_L = (
                - (5.0/48.0)  * cy * w_j * Ec(i,j-3,k) 
                + (13.0/16.0) * cy * w_j * Ec(i,j-2,k) 
                - (17.0/24.0) * cy * w_j * Ec(i,j-1,k)
                - (17.0/24.0) * cy * w_j * Ec(i,j,k)  
                + (13.0/16.0) * cy * w_j * Ec(i,j+1,k)
                - (5.0/48.0)  * cy * w_j * Ec(i,j+2,k))
                / (dy*dy);
            double F_eq_y_L = f_eq_y_L - (dy*dy / 24.0) * d2f_y_L;
            double Sy = (F_eq_y_R - F_eq_y_L) / dy;

            // Z direction
            double f_eq_z_R = cz * w_j * Ef(i,j,k,2);
            double d2f_z_R = (
                - (5.0/48.0)  * cz * w_j * Ec(i,j,k-2) 
                + (13.0/16.0) * cz * w_j * Ec(i,j,k-1) 
                - (17.0/24.0) * cz * w_j * Ec(i,j,k)
                - (17.0/24.0) * cz * w_j * Ec(i,j,k+1) 
                + (13.0/16.0) * cz * w_j * Ec(i,j,k+2)
                - (5.0/48.0)  * cz * w_j * Ec(i,j,k+3))
                / (dz*dz);
            double F_eq_z_R = f_eq_z_R - (dz*dz / 24.0) * d2f_z_R;

            double f_eq_z_L = cz * w_j * Ef(i,j,k-1,2);
            double d2f_z_L = (
                - (5.0/48.0)  * cz * w_j * Ec(i,j,k-3) 
                + (13.0/16.0) * cz * w_j * Ec(i,j,k-2)
                - (17.0/24.0) * cz * w_j * Ec(i,j,k-1)
                - (17.0/24.0) * cz * w_j * Ec(i,j,k)  
                + (13.0/16.0) * cz * w_j * Ec(i,j,k+1)
                - (5.0/48.0)  * cz * w_j * Ec(i,j,k+2))
                / (dz*dz);
            double F_eq_z_L = f_eq_z_L - (dz*dz / 24.0) * d2f_z_L;
            double Sz = (F_eq_z_R - F_eq_z_L) / dz;

            rhs(i,j,k) = -(dFdx_act + dFdy_act + dFdz_act) + (Sx + Sy + Sz);
        });
    }
}

// =========================================================
// 4. Verification Test
// =========================================================
int main(int argc, char* argv[]) {
    amrex::Initialize(argc, argv);
    { 
        amrex::Print() << "Initializing AMReX 3D Literal Local Equilibrium Setup...\n";
        
        amrex::Box domain(amrex::IntVect(0,0,0), amrex::IntVect(N-1,N-1,N-1));
        amrex::RealBox real_box({AMREX_D_DECL(-1.0, -1.0, -1.0)}, {AMREX_D_DECL(1.0, 1.0, 1.0)});
        amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};
        amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

        amrex::BoxArray ba(domain);
        ba.maxSize(8); 
        amrex::DistributionMapping dm(ba);

        // Allocate globals. 3 components (X,Y,Z), 3 ghost cells
        E_center = new amrex::MultiFab(ba, dm, 1, 3);
        E_face   = new amrex::MultiFab(ba, dm, 3, 3); 
        E_prime  = new amrex::MultiFab(ba, dm, 3, 3);

        initialize_equilibrium_arrays(geom);

        // Setup the physical state u
        amrex::MultiFab u(ba, dm, 1, 3);
        amrex::MultiFab::Copy(u, *E_center, 0, 0, 1, 0); 
        u.FillBoundary(geom.periodicity()); 

        // Evaluate RHS
        amrex::MultiFab rhs(ba, dm, 1, 0); 
        get_rhs_local_equilibrium_3D(u, rhs, geom);

        // Track max error
        double max_err = 0.0;
        amrex::Box interior_box(amrex::IntVect(3,3,3), amrex::IntVect(N-4,N-4,N-4));

        for (amrex::MFIter mfi(rhs); mfi.isValid(); ++mfi) {
            amrex::Box bx = mfi.tilebox() & interior_box; 
            if (bx.ok()) {
                auto const& r_arr = rhs.const_array(mfi);
                for(int k=bx.smallEnd(2); k<=bx.bigEnd(2); ++k)
                for(int j=bx.smallEnd(1); j<=bx.bigEnd(1); ++j)
                for(int i=bx.smallEnd(0); i<=bx.bigEnd(0); ++i)
                    max_err = std::max(max_err, std::abs(r_arr(i,j,k)));
            }
        }
        
        amrex::ParallelDescriptor::ReduceRealMax(max_err);

        amrex::Print() << std::scientific << std::setprecision(10);
        amrex::Print() << "--- 3D EXACT STEADY-STATE TEST ---\n"
                       << "Approach: Local Equilibrium Evaluation (AMReX ncomp=3 Parallel)\n"
                       << "Max 3D RHS Error: " << max_err << "\n";
        
        if (max_err < 1e-12) {
            amrex::Print() << "SUCCESS: The scheme is well-balanced to machine precision!\n";
        }

        // Cleanup
        delete E_center; delete E_face; delete E_prime;
    }
    amrex::Finalize();
    return 0;
}
