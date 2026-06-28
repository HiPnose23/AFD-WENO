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
#include <AMReX_Reduce.H> 

#include "../weno_lib.hpp"

// =========================================================
// 1. Setup and Parameters (Isothermal Equilibrium Test)
// =========================================================
constexpr int N = 64; 
constexpr int m = 1; // 1: Cylindrical, 2: Spherical
constexpr double dr = 2.0 / N; 
constexpr double gamma_const = 5.0 / 3.0; 

inline AMREX_GPU_DEVICE double get_r(int i) { return 0.0 + (i + 0.5) * dr; } 

// Arbitrary density profile to show balance holds regardless of mass distribution
inline AMREX_GPU_DEVICE double rho_0(double r) {
    return 1.0 + 0.5 * std::sin(M_PI * r); 
}

// Exact solution is a stationary state (u=0) with uniform pressure (p=1)
inline AMREX_GPU_DEVICE void exact_sol(double r, double t, double& rho, double& rhou, double& E) {
    rho = rho_0(r);
    rhou = 0.0; // Fluid at rest
    double p = 1.0; // Uniform pressure balances the expanding geometry
    E = p / (gamma_const - 1.0); // Kinetic energy is zero
}

// =========================================================
// 2. Boundary Conditions
// =========================================================
void apply_boundaries(amrex::MultiFab& u, const amrex::Geometry& geom) {
    u.FillBoundary(geom.periodicity()); 
    const amrex::Box& domain = geom.Domain();
    int dom_lo = domain.smallEnd(0);
    int dom_hi = domain.bigEnd(0);

#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
    for (amrex::MFIter mfi(u, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const& arr = u.array(mfi);
        amrex::ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            if (i < dom_lo) {
                // Axisymmetric inner boundary (pole)
                int dist = dom_lo - i; 
                int sym_i = dom_lo + dist - 1; 
                arr(i,j,k,0) =  arr(sym_i,j,k,0); 
                arr(i,j,k,1) = -arr(sym_i,j,k,1); // Momentum is odd
                arr(i,j,k,2) =  arr(sym_i,j,k,2); 
            } 
            else if (i > dom_hi) {
                // For the well-balanced test, enforce the exact stationary state
                double r_ghost = get_r(i);
                double rho_e, rhou_e, E_e;
                exact_sol(r_ghost, 0.0, rho_e, rhou_e, E_e);
                arr(i,j,k,0) = rho_e;
                arr(i,j,k,1) = rhou_e;
                arr(i,j,k,2) = E_e;
            }
        });
    }
}

// =========================================================
// 3. AFD-WENO Flux Function (Primitive Reconstruction)
// =========================================================
void compute_flux(const amrex::MultiFab& u, amrex::MultiFab& flux, const amrex::Geometry& geom) {
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
    for (amrex::MFIter mfi(u, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box bx = mfi.tilebox(); 
        amrex::IntVect lo = bx.smallEnd();
        amrex::IntVect hi = bx.bigEnd();
        lo[0] -= 1; 
        amrex::Box flux_box(lo, hi);
        
        auto const& u_arr = u.const_array(mfi);
        auto const& flx = flux.array(mfi);

        amrex::ParallelFor(flux_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            if (i == -1) {
                for (int c = 0; c < 3; ++c) flx(i,j,k,c) = 0.0;
                return;
            }

            auto get_prim = [&](int idx, int comp) {
                double r = u_arr(idx,j,k,0);
                double rhou = u_arr(idx,j,k,1);
                double E = u_arr(idx,j,k,2);
                
                double v = (r > 1e-14) ? (rhou / r) : 0.0;
                double p = (gamma_const - 1.0) * (E - 0.5 * r * v * v);
                
                if (comp == 0) return r;
                if (comp == 1) return v;
                return p;
            };

            double prim_L[3], prim_R[3];

            for (int c = 0; c < 3; ++c) {
                double w_L_i, w_R_i, dwdr_i;
                weno_ao_3_interpolation(
                    get_prim(i-2, c), get_prim(i-1, c), get_prim(i, c), 
                    get_prim(i+1, c), get_prim(i+2, c), dr, w_L_i, w_R_i, dwdr_i
                );
                prim_L[c] = w_R_i;

                double w_L_ip1, w_R_ip1, dwdr_ip1;
                weno_ao_3_interpolation(
                    get_prim(i-1, c), get_prim(i, c), get_prim(i+1, c), 
                    get_prim(i+2, c), get_prim(i+3, c), dr, w_L_ip1, w_R_ip1, dwdr_ip1
                );
                prim_R[c] = w_L_ip1;
            }
            
            auto prim_to_cons = [&](const double* prim, double* cons) {
                double r = prim[0];
                double v = prim[1];
                double p = prim[2];
                cons[0] = r;
                cons[1] = r * v;
                cons[2] = p / (gamma_const - 1.0) + 0.5 * r * v * v;
            };

            double u_L[3], u_R[3];
            prim_to_cons(prim_L, u_L);
            prim_to_cons(prim_R, u_R);

            auto get_phys_flux = [&](const double* state, double* h_out, double& p_out, double& v_out) {
                double rho = state[0], m_mom = state[1], E = state[2];
                v_out = (rho > 1e-14) ? (m_mom / rho) : 0.0;
                p_out = (gamma_const - 1.0) * (E - 0.5 * m_mom * v_out);
                h_out[0] = m_mom;
                h_out[1] = m_mom * v_out + p_out;
                h_out[2] = (E + p_out) * v_out;
            };

            double h_L[3], p_L, v_L; get_phys_flux(u_L, h_L, p_L, v_L);
            double h_R[3], p_R, v_R; get_phys_flux(u_R, h_R, p_R, v_R);

            double c_L = std::sqrt(std::max(0.0, gamma_const * p_L / std::max(u_L[0], 1e-14)));
            double c_R = std::sqrt(std::max(0.0, gamma_const * p_R / std::max(u_R[0], 1e-14)));
            double alpha = std::max(std::abs(v_L) + c_L, std::abs(v_R) + c_R);

            double h_face[3];
            for (int c = 0; c < 3; ++c) {
                h_face[c] = 0.5 * (h_L[c] + h_R[c] - alpha * (u_R[c] - u_L[c]));
            }

            auto F_point = [&](int idx, int comp) {
                double state[3] = {u_arr(idx,j,k,0), u_arr(idx,j,k,1), u_arr(idx,j,k,2)};
                double h_c[3], p_c, v_c; 
                get_phys_flux(state, h_c, p_c, v_c);
                return std::pow(get_r(idx), m) * h_c[comp];
            };
            
            double r_face = get_r(i) + 0.5 * dr;
            double area_face = std::pow(r_face, m);

            for (int c = 0; c < 3; ++c) {
                double d2f = (
                    - (5.0/48.0)  * F_point(i-2, c)
                    + (13.0/16.0) * F_point(i-1, c)
                    - (17.0/24.0) * F_point(i, c)
                    - (17.0/24.0) * F_point(i+1, c)
                    + (13.0/16.0) * F_point(i+2, c)
                    - (5.0/48.0)  * F_point(i+3, c)
                ) / (dr*dr);
                
                flx(i,j,k,c) = (area_face * h_face[c]) - (dr*dr/24.0) * d2f;
            }
        });
    }
}

// =========================================================
// 4. RHS Assembly (With Inline Geometric Divergence Stencil)
// =========================================================
void get_rhs(const amrex::MultiFab& u, amrex::MultiFab& Rhs, const amrex::Geometry& geom) {
    amrex::MultiFab flux_act(u.boxArray(), u.DistributionMap(), 3, 1);
    compute_flux(u, flux_act, geom);

#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
    for (amrex::MFIter mfi(u, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto const& f_act = flux_act.const_array(mfi); 
        auto const& u_arr = u.const_array(mfi);
        auto const& rhs   = Rhs.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double xi = get_r(i); 
            
            double rho = u_arr(i,j,k,0);
            double rhou = u_arr(i,j,k,1);
            double E = u_arr(i,j,k,2);
            double v = (rho > 1e-14) ? (rhou / rho) : 0.0;
            double p = (gamma_const - 1.0) * (E - 0.5 * rho * v * v);

            auto get_f_ref = [&](int idx) {
                double r_f = get_r(idx) + 0.5 * dr;
                double area_f = std::pow(r_f, m);
                double d2f = (
                    - (5.0/48.0)  * std::pow(get_r(idx-2), m)
                    + (13.0/16.0) * std::pow(get_r(idx-1), m)
                    - (17.0/24.0) * std::pow(get_r(idx),   m)
                    - (17.0/24.0) * std::pow(get_r(idx+1), m)
                    + (13.0/16.0) * std::pow(get_r(idx+2), m)
                    - (5.0/48.0)  * std::pow(get_r(idx+3), m)
                ) / (dr*dr);
                return area_f - (dr*dr/24.0) * d2f;
            };

            for (int c = 0; c < 3; ++c) {
                double div_F_act = (f_act(i,j,k,c) - f_act(i-1,j,k,c)) / (std::pow(xi, m) * dr);     
                
                double S = 0.0;
                if (c == 1) {
                    double div_F_ref = (get_f_ref(i) - get_f_ref(i-1)) / (std::pow(xi, m) * dr);
                    S = p * div_F_ref; 
                }
                
                rhs(i,j,k,c) = -div_F_act + S;
            }
        });
    }
}

// =========================================================
// 5. Main Simulation
// =========================================================
int main(int argc, char* argv[]) {
    amrex::Initialize(argc, argv);
    { 
        amrex::Print() << "Initializing Well-Balanced Benchmark (m=" << m << ")...\n";
        
        amrex::Box domain(amrex::IntVect(AMREX_D_DECL(0,0,0)), amrex::IntVect(AMREX_D_DECL(N-1,0,0)));
        amrex::RealBox real_box({AMREX_D_DECL(0.0, -1.0, -1.0)}, {AMREX_D_DECL(2.0, 1.0, 1.0)});
        
        amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
        amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

        amrex::BoxArray ba(domain);
        ba.maxSize(16); 
        amrex::DistributionMapping dm(ba);

        amrex::MultiFab u(ba, dm, 3, 3);
        amrex::MultiFab u1(ba, dm, 3, 3);
        amrex::MultiFab u2(ba, dm, 3, 3);
        amrex::MultiFab rhs(ba, dm, 3, 0);

#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
        for (amrex::MFIter mfi(u, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            auto const& arr = u.array(mfi);
            amrex::ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                double rho_e, rhou_e, E_e;
                exact_sol(get_r(i), 0.0, rho_e, rhou_e, E_e);
                arr(i,j,k,0) = rho_e;
                arr(i,j,k,1) = rhou_e;
                arr(i,j,k,2) = E_e;
            });
        }
        apply_boundaries(u, geom); 

        double t = 0.0, t_end = 0.4, CFL = 0.4;
        int step = 0;

        amrex::Print() << "Starting Time Integration...\n";
        while (t < t_end) {
            // Speed of sound c = sqrt(gamma * p / rho).
            // Here min rho is 0.5, p is 1.0, so max c is approx sqrt(5/3 * 1 / 0.5)
            double max_c = std::sqrt(gamma_const * 1.0 / 0.5); 
            double dt = CFL * dr / max_c; 
            if (t + dt > t_end) dt = t_end - t;

            get_rhs(u, rhs, geom);
            amrex::MultiFab::Copy(u1, u, 0, 0, 3, 0);                 
            amrex::MultiFab::Saxpy(u1, dt, rhs, 0, 0, 3, 0);          
            apply_boundaries(u1, geom);

            get_rhs(u1, rhs, geom);
            amrex::MultiFab::LinComb(u2, 0.75, u, 0, 0.25, u1, 0, 0, 3, 0); 
            amrex::MultiFab::Saxpy(u2, 0.25 * dt, rhs, 0, 0, 3, 0);
            apply_boundaries(u2, geom);

            get_rhs(u2, rhs, geom);
            amrex::MultiFab::LinComb(u, 1.0/3.0, u, 0, 2.0/3.0, u2, 0, 0, 3, 0);
            amrex::MultiFab::Saxpy(u, (2.0/3.0) * dt, rhs, 0, 0, 3, 0);
            apply_boundaries(u, geom);

            t += dt; step++;
            if (step % 10 == 0) amrex::Print() << "Step " << step << " | Time: " << t << " / " << t_end << "\n";
        }

        amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpSum, amrex::ReduceOpSum> reduce_ops;
        amrex::ReduceData<double, double, double> reduce_data(reduce_ops);
        using ReduceTuple = typename amrex::ReduceData<double, double, double>::Type;

#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
        for (amrex::MFIter mfi(u, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            auto const& arr = u.const_array(mfi);
            reduce_ops.eval(mfi.tilebox(), reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
                // For well-balanced property, we check if momentum strictly remained zero
                double err = std::abs(arr(i,j,k,1) - 0.0); 
                
                double r_lo = i * dr, r_hi = (i + 1) * dr;    
                double vol  = (std::pow(r_hi, m + 1.0) - std::pow(r_lo, m + 1.0)) / (m + 1.0);
                return {err, err * vol, err * err * vol};
            });
        }

        auto [linf_err, l1_err, l2_err] = reduce_data.value();
        amrex::ParallelDescriptor::ReduceRealMax(linf_err);
        amrex::ParallelDescriptor::ReduceRealSum(l1_err);
        amrex::ParallelDescriptor::ReduceRealSum(l2_err);
        l2_err = std::sqrt(l2_err);

        amrex::Print() << "\n--- Well-Balanced Verification (Momentum Error) ---\n" << std::scientific << std::setprecision(8)
                       << "L1 Error:    " << l1_err << "\nL2 Error:    " << l2_err << "\nL-inf Error: " << linf_err << "\n";
        
        amrex::WriteSingleLevelPlotfile("plt_wb_final", u, {"rho", "rhou", "E"}, geom, t, step);
    }
    amrex::Finalize();
    return 0;
}
