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

// ==============================================================================
// CONFIG
// ==============================================================================
constexpr int N_R = 256;
constexpr int N_Z = 512;
constexpr double Gamma = 5.0 / 3.0;

// GRID
const double R_min = 0.0, R_max = 10.0;
const double z_min = -10.0, z_max = 10.0;
const double dR = (R_max - R_min) / N_R;
const double dz = (z_max - z_min) / N_Z;

inline AMREX_GPU_DEVICE double get_R(int i) { return R_min + (i + 0.5) * dR; }
inline AMREX_GPU_DEVICE double get_z(int j) { return z_min + (j + 0.5) * dz; }

// ==============================================================================
// CONVERT CONSERVATIVE TO PRIMITIVE VARIABLES AND VICE VERSA
// ==============================================================================
inline AMREX_GPU_DEVICE void cons_to_prim(const double U[4], double V[4]) {
    V[0] = U[0];                 // rho
    V[1] = U[1] / V[0];          // vR
    V[2] = U[2] / V[0];          // vz
    double kin_energy = 0.5 * V[0] * (V[1]*V[1] + V[2]*V[2]);
    V[3] = (Gamma - 1.0) * (U[3] - kin_energy); // p
}

inline AMREX_GPU_DEVICE void prim_to_cons(const double V[4], double U[4]) {
    U[0] = V[0]; // rho
    U[1] = V[0] * V[1]; // rho vR
    U[2] = V[0] * V[2]; // rho vz
    double kin_energy = 0.5 * V[0] * (V[1]*V[1] + V[2]*V[2]);
    U[3] = V[3] / (Gamma - 1.0) + kin_energy; // E
}

// ==============================================================================
// CFL FUNCTION (details in obsidian)
// ==============================================================================
double compute_dt(const amrex::MultiFab& u) {
    amrex::MultiFab speed(u.boxArray(), u.DistributionMap(), 1, 0);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& u_arr = u.const_array(mfi);
        auto const& s_arr = speed.array(mfi);

        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double U_loc[4] = {
                u_arr(i,j,k,0),
                u_arr(i,j,k,1),
                u_arr(i,j,k,2),
                u_arr(i,j,k,3)
            };

            double V_loc[4];
            cons_to_prim(U_loc, V_loc);

            double cs = std::sqrt(Gamma * V_loc[3] / V_loc[0]);
            double max_char = std::max(std::abs(V_loc[1]) + cs,
                                       std::abs(V_loc[2]) + cs);

            s_arr(i,j,k,0) = max_char;
        });
    }

    double max_speed = speed.max(0, 0, false);
    double cfl = 0.4;

    return cfl * std::min(dR, dz) / max_speed;
}

inline AMREX_GPU_DEVICE void calc_flux_R(const double U[4], const double V[4], double F[4]) {
    F[0] = U[0] * V[1];
    F[1] = U[1] * V[1] + V[3];
    F[2] = U[2] * V[1];
    F[3] = (U[3] + V[3]) * V[1];
}

inline AMREX_GPU_DEVICE void calc_flux_Z(const double U[4], const double V[4], double F[4]) {
    F[0] = U[0] * V[2];
    F[1] = U[1] * V[2];
    F[2] = U[2] * V[2] + V[3];
    F[3] = (U[3] + V[3]) * V[2];
}

// ==============================================================================
// INITIAL CONDITIONS & BOUNDARIES (Mignone Spherical Wind)
// ==============================================================================
inline AMREX_GPU_DEVICE void exact_wind_state(double R, double z, double U[4]) {
    double r = std::sqrt(R*R + z*z);
    double V[4];

    if (r <= 1.0) {
        constexpr double c_sw = 3.0e-2;
        const double rr = std::max(r, 1.0e-12);
        const double v_mag = std::tanh(5.0 * rr);

        V[0] = 1.0 / (v_mag * rr * rr);   // rho v r^2 = 1
        V[1] = v_mag * (R / rr);
        V[2] = v_mag * (z / rr);
        V[3] = V[0] * (c_sw * c_sw) / Gamma;
    } else {
        double c_sa = 4.0e-3;
        V[0] = 0.25;
        V[1] = 0.0;
        V[2] = 0.0;
        V[3] = V[0] * (c_sa * c_sa) / Gamma;
    }

    prim_to_cons(V, U);
}

void apply_boundaries(amrex::MultiFab& u, const amrex::Geometry& geom) {
    const amrex::Box& domain = geom.Domain();
    const int dom_lo_R = domain.smallEnd(0);
    const int dom_hi_R = domain.bigEnd(0);
    const int dom_lo_z = domain.smallEnd(1);
    const int dom_hi_z = domain.bigEnd(1);

    constexpr double rho_floor = 1.0e-12;
    constexpr double p_floor   = 1.0e-12;

    // (1) wind region + floors, valid cells only
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& arr = u.array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const double R_val = get_R(i);
            const double z_val = get_z(j);

            if (R_val*R_val + z_val*z_val <= 1.0) {
                double U_wind[4];
                exact_wind_state(R_val, z_val, U_wind);
                for (int n = 0; n < 4; ++n) arr(i,j,k,n) = U_wind[n];
                return;
            }

            double U_loc[4] = { arr(i,j,k,0), arr(i,j,k,1), arr(i,j,k,2), arr(i,j,k,3) };
            double V_loc[4];
            cons_to_prim(U_loc, V_loc);
            if (V_loc[0] < rho_floor || V_loc[3] < p_floor) {
                V_loc[0] = std::max(V_loc[0], rho_floor);
                V_loc[3] = std::max(V_loc[3], p_floor);
                prim_to_cons(V_loc, U_loc);
                for (int n = 0; n < 4; ++n) arr(i,j,k,n) = U_loc[n];
            }
        });
    }

    u.FillBoundary(geom.periodicity());

    // (2) physical ghost cells only: axis reflection + outflow
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& arr = u.array(mfi);
        amrex::ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            if (i >= dom_lo_R && i <= dom_hi_R && j >= dom_lo_z && j <= dom_hi_z) return;

            int i_idx = i;
            double sign_vR = 1.0;
            if (i < dom_lo_R) { i_idx = dom_lo_R + (dom_lo_R - i) - 1; sign_vR = -1.0; }
            else if (i > dom_hi_R) { i_idx = dom_hi_R; }

            int j_idx = j;
            if (j < dom_lo_z) j_idx = dom_lo_z;
            else if (j > dom_hi_z) j_idx = dom_hi_z;

            arr(i,j,k,0) = arr(i_idx,j_idx,k,0);
            arr(i,j,k,1) = arr(i_idx,j_idx,k,1) * sign_vR;
            arr(i,j,k,2) = arr(i_idx,j_idx,k,2);
            arr(i,j,k,3) = arr(i_idx,j_idx,k,3);
        });
    }
}
void compute_flux_R_5(const amrex::MultiFab& u, amrex::MultiFab& flux, const amrex::Geometry& geom) {
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        amrex::Box flux_box = mfi.validbox();
        flux_box.growLo(0, 1);

        auto const& u_arr = u.const_array(mfi);
        auto const& flx = flux.array(mfi);

        amrex::ParallelFor(flux_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double V_L[4], V_R[4];

            for (int n = 0; n < 4; ++n) {
                auto prim_f = [&](int idx) {
                    double U_loc[4] = {
                        u_arr(idx,j,k,0),
                        u_arr(idx,j,k,1),
                        u_arr(idx,j,k,2),
                        u_arr(idx,j,k,3)
                    };
                    double V_loc[4];
                    cons_to_prim(U_loc, V_loc);
                    return V_loc[n];
                };

                double dwdr;
                double dummy_L, dummy_R;
                weno_ao_5_3_interpolation(
                    prim_f(i-2), prim_f(i-1), prim_f(i), prim_f(i+1), prim_f(i+2),
                    dR, dummy_L, V_L[n], dwdr
                );
                weno_ao_5_3_interpolation(
                    prim_f(i-1), prim_f(i), prim_f(i+1), prim_f(i+2), prim_f(i+3),
                    dR, V_R[n], dummy_R, dwdr
                );
            }

            constexpr double rho_floor = 1e-12;
            constexpr double p_floor   = 1e-14;

            // In order to avoid NaNs
            V_L[0] = std::max(V_L[0], rho_floor);
            V_R[0] = std::max(V_R[0], rho_floor);
            V_L[3] = std::max(V_L[3], p_floor);
            V_R[3] = std::max(V_R[3], p_floor);

            double U_L[4], U_R[4], F_L[4], F_R[4];
            prim_to_cons(V_L, U_L);
            prim_to_cons(V_R, U_R);
            calc_flux_R(U_L, V_L, F_L);
            calc_flux_R(U_R, V_R, F_R);

            double cs_L = std::sqrt(Gamma * V_L[3] / V_L[0]);
            double cs_R = std::sqrt(Gamma * V_R[3] / V_R[0]);

            double s_L = std::min(V_L[1] - cs_L, V_R[1] - cs_R);
            double s_R = std::max(V_L[1] + cs_L, V_R[1] + cs_R);


            double F_HLL[4];
            if (s_L >= 0.0) {
                for (int n = 0; n < 4; ++n) F_HLL[n] = F_L[n];
            } else if (s_R <= 0.0) {
                for (int n = 0; n < 4; ++n) F_HLL[n] = F_R[n];
            } else {
                for (int n = 0; n < 4; ++n) {
                    F_HLL[n] = (s_R * F_L[n] - s_L * F_R[n] + s_L * s_R * (U_R[n] - U_L[n])) / (s_R - s_L);
                }
            }

            double R_face = get_R(i) + 0.5 * dR;

            for (int n = 0; n < 4; ++n) {
                auto mod_flux = [&](int idx) {
                    double U_loc[4] = {
                        u_arr(idx,j,k,0),
                        u_arr(idx,j,k,1),
                        u_arr(idx,j,k,2),
                        u_arr(idx,j,k,3)
                    };
                    double V_loc[4];
                    cons_to_prim(U_loc, V_loc);
                    double F_loc[4];
                    calc_flux_R(U_loc, V_loc, F_loc);
                    return get_R(idx) * F_loc[n];
                };

                double d2, d4;
                weno_ao_63_boundary(
                    mod_flux(i-2), mod_flux(i-1), mod_flux(i), mod_flux(i+1), mod_flux(i+2), mod_flux(i+3),
                    dR, d2, d4
                );

                flx(i,j,k,n) = R_face * F_HLL[n]
                             - (dR * dR / 24.0) * d2
                             + (7.0 * std::pow(dR, 4) / 5760.0) * d4;
            }
        });
    }
}

void compute_flux_Z_5(const amrex::MultiFab& u, amrex::MultiFab& flux, const amrex::Geometry& geom) {
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        amrex::Box flux_box = mfi.validbox();
        flux_box.growLo(1, 1);

        auto const& u_arr = u.const_array(mfi);
        auto const& flx = flux.array(mfi);

        amrex::ParallelFor(flux_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double V_L[4], V_R[4];

            for (int n = 0; n < 4; ++n) {
                auto prim_f = [&](int idx) {
                    double U_loc[4] = {
                        u_arr(i,idx,k,0),
                        u_arr(i,idx,k,1),
                        u_arr(i,idx,k,2),
                        u_arr(i,idx,k,3)
                    };
                    double V_loc[4];
                    cons_to_prim(U_loc, V_loc);
                    return V_loc[n];
                };

                double dwdz;
                double dummy_L, dummy_R;
                weno_ao_5_3_interpolation(
                    prim_f(j-2), prim_f(j-1), prim_f(j), prim_f(j+1), prim_f(j+2),
                    dz, dummy_L, V_L[n], dwdz
                );
                weno_ao_5_3_interpolation(
                    prim_f(j-1), prim_f(j), prim_f(j+1), prim_f(j+2), prim_f(j+3),
                    dz, V_R[n], dummy_R, dwdz
                );
            }

            constexpr double rho_floor = 1e-12;
            constexpr double p_floor   = 1e-14;

            // In order to avoid NaNs
            V_L[0] = std::max(V_L[0], rho_floor);
            V_R[0] = std::max(V_R[0], rho_floor);
            V_L[3] = std::max(V_L[3], p_floor);
            V_R[3] = std::max(V_R[3], p_floor);

            double U_L[4], U_R[4], F_L[4], F_R[4];
            prim_to_cons(V_L, U_L);
            prim_to_cons(V_R, U_R);
            calc_flux_Z(U_L, V_L, F_L);
            calc_flux_Z(U_R, V_R, F_R);

            double cs_L = std::sqrt(Gamma * V_L[3] / V_L[0]);
            double cs_R = std::sqrt(Gamma * V_R[3] / V_R[0]);

            double s_L = std::min(V_L[2] - cs_L, V_R[2] - cs_R);
            double s_R = std::max(V_L[2] + cs_L, V_R[2] + cs_R);


            double F_HLL[4];
            if (s_L >= 0.0) {
                for (int n = 0; n < 4; ++n) F_HLL[n] = F_L[n];
            } else if (s_R <= 0.0) {
                for (int n = 0; n < 4; ++n) F_HLL[n] = F_R[n];
            } else {
                for (int n = 0; n < 4; ++n) {
                    F_HLL[n] = (s_R * F_L[n] - s_L * F_R[n] + s_L * s_R * (U_R[n] - U_L[n])) / (s_R - s_L);
                }
            }

            for (int n = 0; n < 4; ++n) {
                auto standard_flux = [&](int idx) {
                    double U_loc[4] = {
                        u_arr(i,idx,k,0),
                        u_arr(i,idx,k,1),
                        u_arr(i,idx,k,2),
                        u_arr(i,idx,k,3)
                    };
                    double V_loc[4];
                    cons_to_prim(U_loc, V_loc);
                    double F_loc[4];
                    calc_flux_Z(U_loc, V_loc, F_loc);
                    return F_loc[n];
                };

                double d2, d4;
                weno_ao_63_boundary(
                    standard_flux(j-2), standard_flux(j-1), standard_flux(j), standard_flux(j+1), standard_flux(j+2), standard_flux(j+3),
                    dz, d2, d4
                );

                flx(i,j,k,n) = F_HLL[n]
                             - (dz * dz / 24.0) * d2
                             + (7.0 * std::pow(dz, 4) / 5760.0) * d4;
            }
        });
    }
}

void get_rhs_2d_euler(const amrex::MultiFab& u, amrex::MultiFab& Rhs, const amrex::Geometry& geom) {
    amrex::MultiFab flux_R(u.boxArray(), u.DistributionMap(), 4, 1);
    amrex::MultiFab flux_Z(u.boxArray(), u.DistributionMap(), 4, 1);

    compute_flux_R_5(u, flux_R, geom);
    compute_flux_Z_5(u, flux_Z, geom);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& u_arr = u.const_array(mfi);
        auto const& fR = flux_R.const_array(mfi);
        auto const& fZ = flux_Z.const_array(mfi);
        auto const& rhs = Rhs.array(mfi);

        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double R_val = get_R(i);

            for (int n = 0; n < 4; ++n) {
                double div_R = (fR(i,j,k,n) - fR(i-1,j,k,n)) / (R_val * dR);
                double div_Z = (fZ(i,j,k,n) - fZ(i,j-1,k,n)) / dz;
                rhs(i,j,k,n) = -div_R - div_Z;
            }

            double U_loc[4] = {
                u_arr(i,j,k,0),
                u_arr(i,j,k,1),
                u_arr(i,j,k,2),
                u_arr(i,j,k,3)
            };
            double V_loc[4];
            cons_to_prim(U_loc, V_loc);

            rhs(i,j,k,1) += V_loc[3] / R_val;
        });
    }
}

// ==============================================================================
// MAIN RUN LOOP (SSPRK5)
// ==============================================================================
int main(int argc, char* argv[]) {
    amrex::Initialize(argc, argv);
    amrex::Print() << "Initializing 2D Cylindrical Euler Solver (AFD-WENO)...\n";

    amrex::Box domain(amrex::IntVect(AMREX_D_DECL(0,0,0)),
                      amrex::IntVect(AMREX_D_DECL(N_R-1, N_Z-1, 0)));
    amrex::RealBox real_box({AMREX_D_DECL(R_min, z_min, 0.0)},
                            {AMREX_D_DECL(R_max, z_max, 1.0)});
    amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
    amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

    amrex::BoxArray ba(domain);
    ba.maxSize(64);
    amrex::DistributionMapping dm(ba);

    amrex::MultiFab u  (ba, dm, 4, 3);
    amrex::MultiFab u0 (ba, dm, 4, 3);
    amrex::MultiFab u1 (ba, dm, 4, 3);
    amrex::MultiFab u2 (ba, dm, 4, 3);
    amrex::MultiFab u3 (ba, dm, 4, 3);
    amrex::MultiFab u4 (ba, dm, 4, 3);
    amrex::MultiFab rhs(ba, dm, 4, 0);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& arr = u.array(mfi);
        amrex::ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double U_init[4];
            exact_wind_state(get_R(i), get_z(j), U_init);
            for (int n = 0; n < 4; ++n) arr(i,j,k,n) = U_init[n];
        });
    }
    apply_boundaries(u, geom);

    double t = 0.0, t_end = 20.0;
    double dt = compute_dt(u);
    bool plotted = 0;

    amrex::Print() << "Starting Time Integration... dt = " << dt << "\n";
    int step = 0;

    while (t < t_end) {
        dt = compute_dt(u);
        if (t + dt > t_end) dt = t_end - t;

        amrex::MultiFab::Copy(u0, u, 0, 0, 4, 3);

        get_rhs_2d_euler(u0, rhs, geom);
        amrex::MultiFab::Copy(u1, u0, 0, 0, 4, 0);
        amrex::MultiFab::Saxpy(u1, dt * 0.39175222700392, rhs, 0, 0, 4, 0);
        apply_boundaries(u1, geom);

        get_rhs_2d_euler(u1, rhs, geom);
        amrex::MultiFab::LinComb(u2, 0.44437049406734, u0, 0,
                                     0.55562950593266, u1, 0, 0, 4, 0);
        amrex::MultiFab::Saxpy(u2, dt * 0.36841059262959, rhs, 0, 0, 4, 0);
        apply_boundaries(u2, geom);

        get_rhs_2d_euler(u2, rhs, geom);
        amrex::MultiFab::LinComb(u3, 0.62010185138540, u0, 0,
                                     0.37989814861460, u2, 0, 0, 4, 0);
        amrex::MultiFab::Saxpy(u3, dt * 0.25189177424738, rhs, 0, 0, 4, 0);
        apply_boundaries(u3, geom);

        amrex::MultiFab rhs_u3(ba, dm, 4, 0);
        get_rhs_2d_euler(u3, rhs_u3, geom);
        amrex::MultiFab::LinComb(u4, 0.17807995410773, u0, 0,
                                     0.82192004589227, u3, 0, 0, 4, 0);
        amrex::MultiFab::Saxpy(u4, dt * 0.54497475021237, rhs_u3, 0, 0, 4, 0);
        apply_boundaries(u4, geom);

        get_rhs_2d_euler(u4, rhs, geom);
        u.setVal(0.0);
        amrex::MultiFab::Saxpy(u, 0.00683325884039, u0, 0, 0, 4, 0);
        amrex::MultiFab::Saxpy(u, 0.51723167208978, u2, 0, 0, 4, 0);
        amrex::MultiFab::Saxpy(u, 0.12759831133288, u3, 0, 0, 4, 0);
        amrex::MultiFab::Saxpy(u, 0.34833675773694, u4, 0, 0, 4, 0);
        amrex::MultiFab::Saxpy(u, dt * 0.08460416338212, rhs_u3, 0, 0, 4, 0);
        amrex::MultiFab::Saxpy(u, dt * 0.22600748319395, rhs, 0, 0, 4, 0);
        apply_boundaries(u, geom);

        if (u.contains_nan(0, u.nComp(), 0)) {
            amrex::Print() << "CRASH: NaNs detected in the domain at step " << step << "!\n";
            break;
        }

        t += dt;
        step++;
        if(t>2 && !plotted) {
            amrex::WriteSingleLevelPlotfile("plt_euler_t2", u, {"rho", "rho_vR", "rho_vz", "E"}, geom, t, step);
            plotted = true;
        }

        if (step % 10 == 0) {
            amrex::Print() << "Step " << step << " | Time: " << t << " / " << t_end << "\n";
        }
    }

    amrex::WriteSingleLevelPlotfile("plt_euler", u, {"rho", "rho_vR", "rho_vz", "E"}, geom, t, step);
    amrex::Finalize();
    return 0;
}
