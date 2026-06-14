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
// 1. Setup and Parameters (1D Cylindrical: r)
// =========================================================
constexpr int N = 128;
constexpr double dr = 2.0 / N;
constexpr double gamma_const = 1.4;

inline double get_r(int i) { return 1.0 + (i + 0.5) * dr; } 

amrex::MultiFab* E_center; // ncomp = 3 (rho, rhou, E)
amrex::MultiFab* E_face;   // ncomp = 3 

void initialize_equilibrium_arrays(const amrex::Geometry& geom) {
    for (amrex::MFIter mfi(*E_center); mfi.isValid(); ++mfi) {
        // USE GROWNTILEBOX TO FILL ALL THE WAY THROUGH GHOST CELLS
        const amrex::Box& bx = mfi.growntilebox(); 
        
        auto const& Ec = E_center->array(mfi);
        auto const& Ef = E_face->array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double r = get_r(i);
            double rho = std::exp(-r); // Isothermal Hydrostatic Balance
            double p = rho; 
            
            Ec(i,j,k,0) = rho;
            Ec(i,j,k,1) = 0.0; // Rest momentum
            Ec(i,j,k,2) = p / (gamma_const - 1.0);
            
            double r_f = r + 0.5 * dr;
            double rho_f = std::exp(-r_f);
            
            Ef(i,j,k,0) = rho_f;
            Ef(i,j,k,1) = 0.0;
            Ef(i,j,k,2) = rho_f / (gamma_const - 1.0);
        });
    }
    // No FillBoundary needed for analytical E arrays
}

// =========================================================
// 2. Custom Hydrostatic Boundary Condition
// =========================================================
void apply_hydrostatic_boundaries_1D(amrex::MultiFab& u, const amrex::Geometry& geom) {
    u.FillBoundary(geom.periodicity());

    const amrex::Box& domain = geom.Domain();
    int dom_lo = domain.smallEnd(0);
    int dom_hi = domain.bigEnd(0);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        const amrex::Box& gbx = mfi.growntilebox(); 
        auto const& arr = u.array(mfi);
        auto const& eq = E_center->array(mfi); 

        amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // Overwrite anything outside the physical domain with exact hydrostatic states
            if (i < dom_lo || i > dom_hi) {
                arr(i,j,k,0) = eq(i,j,k,0); // Matches the exponential drop-off exactly
                arr(i,j,k,1) = eq(i,j,k,1); // Enforces 0.0 resting momentum
                arr(i,j,k,2) = eq(i,j,k,2); // Enforces matching pressure
            } 
        });
    }
}


// =========================================================
// 3. A-WENO Flux Function with Lax-Friedrichs Riemann Solver
// =========================================================
void compute_wb_flux_1D(const amrex::MultiFab& u, amrex::MultiFab& flux, const amrex::Geometry& geom) {
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        amrex::Box bx = mfi.validbox(); 
        
        // Expand the flux box 1 cell left to populate f_act(-1) for the divergence
        amrex::IntVect lo = bx.smallEnd();
        amrex::IntVect hi = bx.bigEnd();
        lo[0] -= 1; 
        amrex::Box flux_box(lo, hi);
        
        auto const& u_arr = u.const_array(mfi);
        auto const& Ec = E_center->const_array(mfi);
        auto const& Ef = E_face->const_array(mfi); 
        auto const& flx = flux.array(mfi);

        amrex::ParallelFor(flux_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double u_L[3], u_R[3];

            // 1. Reconstruct Left and Right States at the face (i + 1/2)
            for (int c = 0; c < 3; ++c) {
                double w_L_i, w_R_i, dwdr_i;
                // Left-biased stencil (centered on cell i)
                weno_ao_3_interpolation(
                    u_arr(i-2,j,k,c) - Ec(i-2,j,k,c), u_arr(i-1,j,k,c) - Ec(i-1,j,k,c),
                    u_arr(i,j,k,c)   - Ec(i,j,k,c),   u_arr(i+1,j,k,c) - Ec(i+1,j,k,c), 
                    u_arr(i+2,j,k,c) - Ec(i+2,j,k,c), dr, w_L_i, w_R_i, dwdr_i
                );
                u_L[c] = w_R_i + Ef(i,j,k,c); // w_R is the right face of cell i

                double w_L_ip1, w_R_ip1, dwdr_ip1;
                // Right-biased stencil (centered on cell i+1)
                weno_ao_3_interpolation(
                    u_arr(i-1,j,k,c) - Ec(i-1,j,k,c), u_arr(i,j,k,c)   - Ec(i,j,k,c),
                    u_arr(i+1,j,k,c) - Ec(i+1,j,k,c), u_arr(i+2,j,k,c) - Ec(i+2,j,k,c), 
                    u_arr(i+3,j,k,c) - Ec(i+3,j,k,c), dr, w_L_ip1, w_R_ip1, dwdr_ip1
                );
                u_R[c] = w_L_ip1 + Ef(i,j,k,c); // w_L is the left face of cell i+1
            }
            
            // 2. Compute Physical Fluxes for Left and Right states
            double rho_L = u_L[0], m_L = u_L[1], E_L = u_L[2];
            double u_vel_L = (rho_L > 1e-14) ? (m_L / rho_L) : 0.0;
            double p_L = (gamma_const - 1.0) * (E_L - 0.5 * m_L * u_vel_L);
            double h_L[3] = {m_L, m_L * u_vel_L + p_L, (E_L + p_L) * u_vel_L};

            double rho_R = u_R[0], m_R = u_R[1], E_R = u_R[2];
            double u_vel_R = (rho_R > 1e-14) ? (m_R / rho_R) : 0.0;
            double p_R = (gamma_const - 1.0) * (E_R - 0.5 * m_R * u_vel_R);
            double h_R[3] = {m_R, m_R * u_vel_R + p_R, (E_R + p_R) * u_vel_R};

            // 3. Local Lax-Friedrichs Flux Splitting
            double c_L = std::sqrt(std::max(0.0, gamma_const * p_L / rho_L));
            double c_R = std::sqrt(std::max(0.0, gamma_const * p_R / rho_R));
            double alpha = std::max(std::abs(u_vel_L) + c_L, std::abs(u_vel_R) + c_R);

            double h[3];
            for (int c = 0; c < 3; ++c) {
                h[c] = 0.5 * (h_L[c] + h_R[c] - alpha * (u_R[c] - u_L[c]));
            }

            // Helper lambda for spatial correction
            auto F= [&](int idx, int comp) {
                double rho_c = u_arr(idx,j,k,0), m_c = u_arr(idx,j,k,1), eng_c = u_arr(idx,j,k,2);
                double uv = (rho_c > 1e-14) ? (m_c / rho_c) : 0.0;
                double pr = (gamma_const - 1.0) * (eng_c - 0.5 * m_c * uv);
                if (comp == 0) return m_c;
                if (comp == 1) return m_c * uv + pr;
                return (eng_c + pr) * uv; 
            };
            
            double r_face = get_r(i) + 0.5 * dr;

            // 4. 6-Point Spatial Correction
            for (int c = 0; c < 3; ++c) {
                double d2f = (
                    - (5.0/48.0)  * get_r(i-2) * F(i-2, c)
                    + (13.0/16.0) * get_r(i-1) * F(i-1, c)
                    - (17.0/24.0) * get_r(i)   * F(i, c)
                    - (17.0/24.0) * get_r(i+1) * F(i+1, c)
                    + (13.0/16.0) * get_r(i+2) * F(i+2, c)
                    - (5.0/48.0)  * get_r(i+3) * F(i+3, c)
                ) / (dr*dr);
                
                flx(i,j,k,c) = (r_face * h[c]) - (dr*dr/24.0) * d2f;
            }
        });
    }
}
// =========================================================
// 4. RHS Assembly 
// =========================================================
void get_rhs_local_equilibrium_1D(const amrex::MultiFab& u, amrex::MultiFab& Rhs, const amrex::Geometry& geom) {
    amrex::MultiFab flux_act(u.boxArray(), u.DistributionMap(), 3, 1);
    compute_wb_flux_1D(u, flux_act, geom);

    amrex::MultiFab flux_eq(u.boxArray(), u.DistributionMap(), 3, 1);
    compute_wb_flux_1D(*E_center, flux_eq, geom);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto const& f_act = flux_act.const_array(mfi); 
        auto const& f_eq  = flux_eq.const_array(mfi);
        auto const& u_arr = u.const_array(mfi);
        auto const& Ec    = E_center->const_array(mfi);
        auto const& rhs   = Rhs.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double r_i = get_r(i);
            for (int c = 0; c < 3; ++c) {
                double dFdr_act = (f_act(i,j,k,c) - f_act(i-1,j,k,c)) / (r_i * dr);
                double Div_E_r  = (f_eq(i,j,k,c)  - f_eq(i-1,j,k,c))  / (r_i * dr);

                double Sr = 0.0;
                if (c == 1) { // Momentum scaled strictly by density ratio
                    Sr = (u_arr(i,j,k,0) / Ec(i,j,k,0)) * Div_E_r;
                }
                if (std::isnan(dFdr_act) || std::isnan(Sr) || std::isnan(u_arr(i,j,k,0))) {
                    printf("NaN detected at i=%d, c=%d | u_rho=%.6e, Ec_rho=%.6e, Div_E_r=%.6e\n", 
                            i, c, u_arr(i,j,k,0), Ec(i,j,k,0), Div_E_r);
                    amrex::Abort("NaN");
                }
                rhs(i,j,k,c) = -dFdr_act + Sr;
            }
        });
    }
}

// =========================================================
// 5. Main Simulation & SSP-RK3
// =========================================================
int main(int argc, char* argv[]) {
    amrex::Initialize(argc, argv);
    { 
        amrex::Print() << "Initializing AMReX A-WENO 1D Cylindrical Euler Solver...\n";
        
        amrex::Box domain(amrex::IntVect(AMREX_D_DECL(0,0,0)), amrex::IntVect(AMREX_D_DECL(N-1,0,0)));
        amrex::RealBox real_box({AMREX_D_DECL(1.0, 0.0, 0.0)}, {AMREX_D_DECL(3.0, 1.0, 1.0)});
        
        // 1. Turn OFF domain periodicity so custom hydrostatic edges apply
        amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
        amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

        amrex::BoxArray ba(domain);
        ba.maxSize(128); 
        amrex::DistributionMapping dm(ba);

        E_center = new amrex::MultiFab(ba, dm, 3, 3);
        E_face   = new amrex::MultiFab(ba, dm, 3, 3); 
        initialize_equilibrium_arrays(geom);

        amrex::MultiFab u(ba, dm, 3, 3);
        amrex::MultiFab u1(ba, dm, 3, 3);
        amrex::MultiFab u2(ba, dm, 3, 3);
        amrex::MultiFab rhs(ba, dm, 3, 0);

        amrex::MultiFab::Copy(u, *E_center, 0, 0, 3, 3); 
        apply_hydrostatic_boundaries_1D(u, geom); 

        double t = 0.0, t_end = 1.0, CFL = 0.4;

        int step = 0;
        while (t < t_end) {
            double max_c = 2.0; 
            double dt = CFL * dr / max_c;
            if (t + dt > t_end) dt = t_end - t;

            get_rhs_local_equilibrium_1D(u, rhs, geom);
            amrex::MultiFab::Copy(u1, u, 0, 0, 3, 0);                 
            amrex::MultiFab::Saxpy(u1, dt, rhs, 0, 0, 3, 0);          
            apply_hydrostatic_boundaries_1D(u1, geom);

            get_rhs_local_equilibrium_1D(u1, rhs, geom);
            amrex::MultiFab::LinComb(u2, 0.75, u, 0, 0.25, u1, 0, 0, 3, 0); 
            amrex::MultiFab::Saxpy(u2, 0.25 * dt, rhs, 0, 0, 3, 0);
            apply_hydrostatic_boundaries_1D(u2, geom);

            get_rhs_local_equilibrium_1D(u2, rhs, geom);
            amrex::MultiFab::LinComb(u, 1.0/3.0, u, 0, 2.0/3.0, u2, 0, 0, 3, 0);
            amrex::MultiFab::Saxpy(u, (2.0/3.0) * dt, rhs, 0, 0, 3, 0);
            apply_hydrostatic_boundaries_1D(u, geom);

            t += dt;
            step++;
            if (step % 10 == 0) amrex::Print() << "Step " << step << " | Time: " << t << " / " << t_end << "\n";
        }

        amrex::MultiFab error(ba, dm, 3, 0);
        amrex::MultiFab::Copy(error, u, 0, 0, 3, 0);
        amrex::MultiFab::Subtract(error, *E_center, 0, 0, 3, 0);

        amrex::Print() << "\n--- Well-Balanced Verification (Density) ---\n"
                       << std::scientific << std::setprecision(8)
                       << "L1 Error:    " << error.norm1(0) * dr << "\n"
                       << "L-inf Error: " << error.norm0(0) << "\n"
                       << "----------------------------------\n";
        
        amrex::WriteSingleLevelPlotfile("plt_final", u, {"rho", "rhou", "E"}, geom, t, step);
        delete E_center; delete E_face;
    }
    amrex::Finalize();
    return 0;
}
