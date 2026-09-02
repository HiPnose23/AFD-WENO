// Note:
// Balance law solved here:
//
//   CARTESIAN   d_t u + d_z( alpha_v u )                 = S(u,z)
//   RADIAL      d_t u + (1/r^m) d_r( r^(m+1) alpha_v u ) = S(u,r)
//
// with the source chosen so that a prescribed profile is an exact steady state:
//
//   CARTESIAN   S(u,z) = -k_eq u                        =>  u_e(z) = exp(-k_eq z / alpha_v)
//   RADIAL      S(u,r) = alpha_v[(m+1) - 2 k_eq r^2] u  =>  u_e(r) = exp(-k_eq r^2)
//
// In both cases S is proportional to u, so S(u)/S(u_e) = u/u_e.
//
// =============================================================================
// How to use
// =============================================================================
// Set mode to either Cartesian or Radial
// Set m if mode is RADIAL. 1: Cylindrical case, 2: Spherical case
// PERT_AMP is the amplitutde of the gaussian perturbation, PERT_X is where it is centered
// This code will produce the results for the non-WB AND the WB scheme.
// ==============================================================================

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>

#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Array.H>
#include <AMReX_RealBox.H>

#include "../weno_lib.hpp"

// ==============================================================================
// CONFIG
// ==============================================================================
enum class Mode { CARTESIAN, RADIAL };
constexpr Mode TEST_MODE = Mode::RADIAL;

constexpr int m = 2; // 1 for Cylindrical, 2 for Spherical (RADIAL only)

constexpr int N = 128;
constexpr double alpha_v = 1.0;

constexpr int p_r = (m == 2) ? 1 : (m + 1);

constexpr double k_eq = 5.0;

constexpr double PERT_AMP = 1.0e-9;
constexpr double PERT_X = 0.5;

constexpr double T_END = 0.3;

// GRID
const double dx = 2.0 / N;

inline AMREX_GPU_DEVICE double get_x(int i) { return 0.0 + (i + 0.5) * dx; }

// ==============================================================================
// EQUILIBRIUM AND SOURCE
// ==============================================================================
inline AMREX_GPU_DEVICE double u_equilibrium(double x) {
    if constexpr (TEST_MODE == Mode::CARTESIAN) {
        return std::exp(-k_eq * x / alpha_v);
    } else {
        return std::exp(-k_eq * x * x);
    }
}

inline AMREX_GPU_DEVICE double source_S(double u, double x) {
    if constexpr (TEST_MODE == Mode::CARTESIAN) {
        return -k_eq * u;
    } else {
        return alpha_v * ((m + 1.0) - 2.0 * k_eq * x * x) * u;
    }
}

inline AMREX_GPU_DEVICE double source_ratio(double u, double u_e) {
    return u / u_e;
}

inline AMREX_GPU_DEVICE double u_initial(double x) {
    return u_equilibrium(x) + PERT_AMP * std::exp(-100.0 * (x - PERT_X) * (x - PERT_X));
}

// ==============================================================================
// BOUNDARIES
// ==============================================================================

void apply_boundaries(amrex::MultiFab& u, const amrex::Geometry& geom) {
    u.FillBoundary(geom.periodicity());

    const amrex::Box& domain = geom.Domain();
    int dom_lo = domain.smallEnd(0);
    int dom_hi = domain.bigEnd(0);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& arr = u.array(mfi);
        amrex::Box gbx = mfi.validbox();
        gbx.grow(0, u.nGrow());

        amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            if (i < dom_lo) {
                if constexpr (TEST_MODE == Mode::RADIAL) {
                    int dist = dom_lo - i;
                    arr(i,j,k) = arr(dom_lo + dist - 1, j, k);
                } else {
                    arr(i,j,k) = u_equilibrium(get_x(i));
                }
            } else if (i > dom_hi) {
                arr(i,j,k) = u_equilibrium(get_x(i));
            }
        });
    }
}


void compute_flux_5(const amrex::MultiFab& u, amrex::MultiFab& flux, const amrex::Geometry& geom) {
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        amrex::Box flux_box = mfi.validbox();
        flux_box.growLo(0, 1);

        auto const& u_arr = u.const_array(mfi);
        auto const& flx = flux.array(mfi);

        amrex::ParallelFor(flux_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double w_L, w_R, dwdx;
            auto point_f = [&](int idx) {
                if constexpr (TEST_MODE == Mode::CARTESIAN) {
                    return alpha_v * u_arr(idx,j,k);
                } else {
                    return std::pow(get_x(idx), p_r) * alpha_v * u_arr(idx,j,k);
                }
            };

            weno_ao_5_3_interpolation(u_arr(i-2,j,k), u_arr(i-1,j,k), u_arr(i,j,k), u_arr(i+1,j,k), u_arr(i+2,j,k), dx, w_L, w_R, dwdx);

            double d2, d4;
            weno_ao_63_boundary(point_f(i-2), point_f(i-1), point_f(i), point_f(i+1), point_f(i+2), point_f(i+3), dx, d2, d4);

            double corr = -(dx * dx / 24.0) * d2 + (7.0 * dx * dx * dx * dx / 5760.0) * d4;

            if constexpr (TEST_MODE == Mode::CARTESIAN) {
                flx(i,j,k,0) = alpha_v * w_R + corr;
            } else {
                double x_face = get_x(i) + 0.5 * dx;
                flx(i,j,k,0) = std::pow(x_face, p_r) * alpha_v * w_R + corr;
            }
        });
    }
}
// D from Section 3.5 in Thesis paper
void apply_operator_5(const amrex::MultiFab& u, amrex::MultiFab& D, const amrex::Geometry& geom) {
    amrex::MultiFab flux_act(u.boxArray(), u.DistributionMap(), 1, 1);
    flux_act.setVal(0.0);
    compute_flux_5(u, flux_act, geom);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& u_arr = u.const_array(mfi);
        auto const& f_act = flux_act.const_array(mfi);
        auto const& d_arr = D.array(mfi);

        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double div_F = (f_act(i,j,k,0) - f_act(i-1,j,k,0)) / dx;

            if constexpr (TEST_MODE == Mode::CARTESIAN) {
                d_arr(i,j,k) = div_F;
            } else if constexpr (m == 2) {
                d_arr(i,j,k) = div_F + m * alpha_v * u_arr(i,j,k);
            } else {
                d_arr(i,j,k) = (f_act(i,j,k,0) - f_act(i-1,j,k,0)) / (std::pow(get_x(i), m) * dx);
            }
        });
    }
}

// ==============================================================================
// RIGHT HAND SIDE
// ==============================================================================
// Well balanced:  d_t u_i = -D[u]_i + ( S(u_i)/S(u_e,i) ) D[u_e]_i
// Standard:       d_t u_i = -D[u]_i + S(u_i, x_i)

void get_rhs_5(const amrex::MultiFab& u, const amrex::MultiFab& u_eq,
               const amrex::MultiFab& D_eq, amrex::MultiFab& Rhs,
               bool use_wb, const amrex::Geometry& geom) {
    amrex::MultiFab D_act(u.boxArray(), u.DistributionMap(), 1, 0);
    D_act.setVal(0.0);
    apply_operator_5(u, D_act, geom);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& u_arr  = u.const_array(mfi);
        auto const& ue_arr = u_eq.const_array(mfi);
        auto const& d_act  = D_act.const_array(mfi);
        auto const& d_eq   = D_eq.const_array(mfi);
        auto const& rhs    = Rhs.array(mfi);

        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            if (use_wb) {
                rhs(i,j,k) = -d_act(i,j,k) + source_ratio(u_arr(i,j,k), ue_arr(i,j,k)) * d_eq(i,j,k);
            } else {
                rhs(i,j,k) = -d_act(i,j,k) + source_S(u_arr(i,j,k), get_x(i));
            }
        });
    }
}

// ==============================================================================
// DEVIATION FROM EQUILIBRIUM
// ==============================================================================

void deviation_norms(const amrex::MultiFab& u, const amrex::BoxArray& ba,
                     const amrex::DistributionMapping& dm,
                     double& linf, double& l1) {
    amrex::MultiFab dev(ba, dm, 2, 0);
    dev.setVal(0.0);

    for (amrex::MFIter mfi(dev); mfi.isValid(); ++mfi) {
        auto const& u_arr = u.const_array(mfi);
        auto const& d_arr = dev.array(mfi);

        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double err = std::abs(u_arr(i,j,k) - u_equilibrium(get_x(i)));
            double vol = dx;
            if constexpr (TEST_MODE == Mode::RADIAL) {
                double x_lo = i * dx;
                double x_hi = (i + 1) * dx;
                vol = (std::pow(x_hi, m + 1.0) - std::pow(x_lo, m + 1.0)) / (m + 1.0);
            }
            d_arr(i,j,k,0) = err;
            d_arr(i,j,k,1) = err * vol;
        });
    }

    linf = dev.norm0(0);
    l1   = dev.sum(1);
}

// ==============================================================================
// RUN
// ==============================================================================

void run_scheme_5(bool use_wb, const amrex::BoxArray& ba,
                  const amrex::DistributionMapping& dm,
                  const amrex::Geometry& geom,
                  const amrex::MultiFab& u_eq, const amrex::MultiFab& D_eq,
                  const std::string& plotname) {
    amrex::MultiFab u (ba, dm, 1, 3);
    amrex::MultiFab u0(ba, dm, 1, 3);
    amrex::MultiFab u1(ba, dm, 1, 3);
    amrex::MultiFab u2(ba, dm, 1, 3);
    amrex::MultiFab u3(ba, dm, 1, 3);
    amrex::MultiFab u4(ba, dm, 1, 3);
    amrex::MultiFab rhs(ba, dm, 1, 0);
    amrex::MultiFab rhs_u3(ba, dm, 1, 0);

    u.setVal(0.0);
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& arr = u.array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            arr(i,j,k) = u_initial(get_x(i));
        });
    }
    apply_boundaries(u, geom);
 
    double t = 0.0;
    double c_speed = (TEST_MODE == Mode::CARTESIAN) ? alpha_v : (alpha_v * 2.0);
    double dt = 0.4 * dx / c_speed;
    int step = 0;

    while (t < T_END) {
        if (t + dt > T_END) dt = T_END - t;
        if (dt <= 0.0) break;

        amrex::MultiFab::Copy(u0, u, 0, 0, 1, 3);

        // SSPRK Stage 1
        get_rhs_5(u0, u_eq, D_eq, rhs, use_wb, geom);
        amrex::MultiFab::Copy(u1, u0, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u1, dt * 0.39175222700392, rhs, 0, 0, 1, 0);
        apply_boundaries(u1, geom);

        // Stage 2
        get_rhs_5(u1, u_eq, D_eq, rhs, use_wb, geom);
        amrex::MultiFab::LinComb(u2, 0.44437049406734, u0, 0, 0.55562950593266, u1, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u2, dt * 0.36841059262959, rhs, 0, 0, 1, 0);
        apply_boundaries(u2, geom);

        // Stage 3
        get_rhs_5(u2, u_eq, D_eq, rhs, use_wb, geom);
        amrex::MultiFab::LinComb(u3, 0.62010185138540, u0, 0, 0.37989814861460, u2, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u3, dt * 0.25189177424738, rhs, 0, 0, 1, 0);
        apply_boundaries(u3, geom);

        // Stage 4
        get_rhs_5(u3, u_eq, D_eq, rhs_u3, use_wb, geom);
        amrex::MultiFab::LinComb(u4, 0.17807995410773, u0, 0, 0.82192004589227, u3, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u4, dt * 0.54497475021237, rhs_u3, 0, 0, 1, 0);
        apply_boundaries(u4, geom);

        // Stage 5
        get_rhs_5(u4, u_eq, D_eq, rhs, use_wb, geom);
        u.setVal(0.0);
        amrex::MultiFab::Saxpy(u, 0.00683325884039, u0, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u, 0.51723167208978, u2, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u, 0.12759831133288, u3, 0, 0, 1, 0);        
        amrex::MultiFab::Saxpy(u, 0.34833675773694, u4, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u, dt * 0.08460416338212, rhs_u3, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u, dt * 0.22600748319395, rhs, 0, 0, 1, 0);
        apply_boundaries(u, geom);

        t += dt; step++;
        if(step==1){
            double linf, l1;
            deviation_norms(u, ba, dm, linf, l1);

            amrex::Print() << std::scientific << std::setprecision(8)
                   << "  steps: " << step << "\n"
                   << "  L1  |u-u_e|: " << l1 << "\n"
                   << "  Linf|u-u_e|: " << linf << "\n";

        }
    }

    double linf, l1;
    deviation_norms(u, ba, dm, linf, l1);

    amrex::Print() << std::scientific << std::setprecision(8)
                   << "  steps: " << step << "\n"
                   << "  L1  |u-u_e|: " << l1 << "\n"
                   << "  Linf|u-u_e|: " << linf << "\n";

    amrex::MultiFab out(ba, dm, 3, 0);
    out.setVal(0.0);
    for (amrex::MFIter mfi(out); mfi.isValid(); ++mfi) {
        auto const& u_arr = u.const_array(mfi);
        auto const& o_arr = out.array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double ue = u_equilibrium(get_x(i));
            o_arr(i,j,k,0) = u_arr(i,j,k);
            o_arr(i,j,k,1) = ue;
            o_arr(i,j,k,2) = u_arr(i,j,k) - ue;
        });
    }
    amrex::WriteSingleLevelPlotfile(plotname, out, {"u", "u_eq", "dev"}, geom, t, step);
}

// ==============================================================================
// MAIN
// ==============================================================================

int main(int argc, char* argv[]) {
    amrex::Initialize(argc, argv);
    {       
        amrex::Print() << "Initializing AMReX Well-Balanced Test (5th Order)...\n";
        amrex::Print() << "Mode: " << ((TEST_MODE == Mode::CARTESIAN) ? "Cartesian (z)" : "Radial") << "\n";
        if constexpr (TEST_MODE == Mode::RADIAL)
            amrex::Print() << "Geometry: " << ((m == 1) ? "Cylindrical (m=1)" : "Spherical (m=2)") << "\n";
        amrex::Print() << "N = " << N << "  k_eq = " << k_eq << "  perturbation = " << PERT_AMP << "\n";

        amrex::Box domain(amrex::IntVect(AMREX_D_DECL(0,0,0)), amrex::IntVect(AMREX_D_DECL(N-1,0,0)));
        amrex::RealBox real_box({AMREX_D_DECL(0.0, -1.0, -1.0)}, {AMREX_D_DECL(2.0, 1.0, 1.0)});
        amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
        amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

        amrex::BoxArray ba(domain);
        ba.maxSize(32);
        amrex::DistributionMapping dm(ba);

        amrex::MultiFab u_eq(ba, dm, 1, 3);
        u_eq.setVal(0.0);
        for (amrex::MFIter mfi(u_eq); mfi.isValid(); ++mfi) {
            auto const& arr = u_eq.array(mfi);
            amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                arr(i,j,k) = u_equilibrium(get_x(i));
            });
        }
        apply_boundaries(u_eq, geom);

        amrex::MultiFab D_eq(ba, dm, 1, 0);
        D_eq.setVal(0.0);
        apply_operator_5(u_eq, D_eq, geom); // compute the operator of equilibrium state only once as it remains unchanged
 
        amrex::MultiFab res(ba, dm, 1, 0);
        get_rhs_5(u_eq, u_eq, D_eq, res, false, geom);
        amrex::Print() << "residual at equilibrium, standard      = " << res.norm0() << "\n";
        get_rhs_5(u_eq, u_eq, D_eq, res, true, geom);
        amrex::Print() << "residual at equilibrium, well-balanced = " << res.norm0() << "\n";

        amrex::Print() << "\n--- Standard scheme ---\n";
        run_scheme_5(false, ba, dm, geom, u_eq, D_eq, "plt_nowb");

        amrex::Print() << "\n--- Well-balanced scheme ---\n";
        run_scheme_5(true, ba, dm, geom, u_eq, D_eq, "plt_wb");
    }
    amrex::Finalize();
    return 0;
}
