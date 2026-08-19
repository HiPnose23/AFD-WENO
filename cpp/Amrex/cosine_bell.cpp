// Note:
// Here we set the axes up as follows:
// i = theta
// j = phi

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
#include <AMReX_Reduce.H>

#include "../weno_lib.hpp"

// ==============================================================================
// CONFIG
// ==============================================================================
constexpr int N_THETA = 32;
constexpr int N_PHI = 2 * N_THETA;

constexpr double u0 = 2.0 * M_PI;
constexpr double sigma_0 = 1.0 / 3.0;
constexpr double theta_c = M_PI / 2.0;
constexpr double phi_c = 3.0 * M_PI / 2.0;

// GRID
const double dtheta = M_PI / N_THETA;
const double dphi = (2.0 * M_PI) / N_PHI;

inline AMREX_GPU_DEVICE double get_theta(int i) { return 0.0 + (i + 0.5) * dtheta; }
inline AMREX_GPU_DEVICE double get_phi(int j) { return 0.0 + (j + 0.5) * dphi; }

// ==============================================================================
// VELOCITY FIELD
// ==============================================================================
inline AMREX_GPU_DEVICE double vel_theta(double theta, double phi) {
    return u0 * std::sin(phi);
}

inline AMREX_GPU_DEVICE double vel_phi(double theta, double phi) {
    return u0 * std::cos(theta) * std::cos(phi);
}

// ==============================================================================
// EXACT SOLUTION (= Initial condition) WILL BE 2ND ORDER AS NOT C^2
// ==============================================================================
inline AMREX_GPU_DEVICE double exact_sol_bell(double theta, double phi) {
    double dot = std::cos(theta_c) * std::cos(theta)
               + std::sin(theta_c) * std::sin(theta) * std::cos(phi - phi_c);
    if (dot > 1.0) dot = 1.0;
    if (dot < -1.0) dot = -1.0;
    double sigma = std::acos(dot);

    if (sigma < sigma_0) return 0.5 * (1.0 + std::cos(M_PI * sigma / sigma_0));
    return 0.0;
}

// Phi is periodic so simple, theta needs to be handled differently.
// Imagine a sphere and you are going up to the north pole on a certain longitutde (phi), lets say 0, you reach the pole and continue going forward,
// you are now going south but on the exact opposite longitude (pi)
void apply_boundaries(amrex::MultiFab& u, const amrex::Geometry& geom) {
    u.FillBoundary(geom.periodicity());

    const amrex::Box& domain = geom.Domain();
    int dom_lo_t = domain.smallEnd(0);
    int dom_hi_t = domain.bigEnd(0);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& arr = u.array(mfi);
        amrex::Box grown_box = mfi.growntilebox();

        amrex::ParallelFor(grown_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            if (i < dom_lo_t || i > dom_hi_t) {
                int i_idx = (i < dom_lo_t) ? (dom_lo_t + (dom_lo_t - i) - 1)
                                           : (dom_hi_t - (i - dom_hi_t) + 1);
                int j_idx = ((j + N_PHI / 2) % N_PHI + N_PHI) % N_PHI; // this formula because if j < 0 we will get negative moduli
                arr(i, j, k) = arr(i_idx, j_idx, k);
            }
        });
    }
}

void compute_flux_theta_5(const amrex::MultiFab& u, amrex::MultiFab& flux, const amrex::Geometry& geom) {
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        amrex::Box flux_box = mfi.validbox();
        flux_box.growLo(0, 1);

        auto const& u_arr = u.const_array(mfi);
        auto const& flx = flux.array(mfi);

        amrex::ParallelFor(flux_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double w_L, w_R, dw, dummy_L, dummy_R;
            auto point_f = [&](int idx) {
                double th = get_theta(idx);
                return std::sin(th) * vel_theta(th, get_phi(j)) * u_arr(idx,j,k);
            };

            weno_ao_5_3_interpolation(u_arr(i-2,j,k), u_arr(i-1,j,k), u_arr(i,j,k), u_arr(i+1,j,k), u_arr(i+2,j,k), dtheta, dummy_L, w_R, dw);
            weno_ao_5_3_interpolation(u_arr(i-1,j,k), u_arr(i,j,k), u_arr(i+1,j,k), u_arr(i+2,j,k), u_arr(i+3,j,k), dtheta, w_L, dummy_R, dw);

            double d2, d4;
            weno_ao_63_boundary(point_f(i-2), point_f(i-1), point_f(i), point_f(i+1), point_f(i+2), point_f(i+3), dtheta, d2, d4);

            double t_face = get_theta(i) + 0.5 * dtheta;
            double v_face = vel_theta(t_face, get_phi(j));
            double q_face = (v_face > 0.0) ? w_R : w_L;

            flx(i,j,k,0) = std::sin(t_face) * v_face * q_face - (dtheta * dtheta / 24.0) * d2 + (7.0 * std::pow(dtheta,4) / 5760.0) * d4;
        });
    }
}

void compute_flux_phi_5(const amrex::MultiFab& u, amrex::MultiFab& flux, const amrex::Geometry& geom) {
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        amrex::Box flux_box = mfi.validbox();
        flux_box.growLo(1, 1);

        auto const& u_arr = u.const_array(mfi);
        auto const& flx = flux.array(mfi);

        amrex::ParallelFor(flux_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double w_L, w_R, dw, dummy_L, dummy_R;
            auto point_f = [&](int idx) {
                double th = get_theta(i);
                return vel_phi(th, get_phi(idx)) * u_arr(i,idx,k);
            };

            weno_ao_5_3_interpolation(u_arr(i,j-2,k), u_arr(i,j-1,k), u_arr(i,j,k), u_arr(i,j+1,k), u_arr(i,j+2,k), dphi, dummy_L, w_R, dw);
            weno_ao_5_3_interpolation(u_arr(i,j-1,k), u_arr(i,j,k), u_arr(i,j+1,k), u_arr(i,j+2,k), u_arr(i,j+3,k), dphi, w_L, dummy_R, dw);

            double d2, d4;
            weno_ao_63_boundary(point_f(j-2), point_f(j-1), point_f(j), point_f(j+1), point_f(j+2), point_f(j+3), dphi, d2, d4);

            double p_face = get_phi(j) + 0.5 * dphi;
            double v_face = vel_phi(get_theta(i), p_face);
            double q_face = (v_face > 0.0) ? w_R : w_L;

            flx(i,j,k,0) = v_face * q_face - (dphi * dphi / 24.0) * d2 + (7.0 * std::pow(dphi,4) / 5760.0) * d4;
        });
    }
}

void get_rhs_bell_5(const amrex::MultiFab& u, amrex::MultiFab& Rhs, const amrex::Geometry& geom) {
    amrex::MultiFab flux_t(u.boxArray(), u.DistributionMap(), 1, 1);
    amrex::MultiFab flux_p(u.boxArray(), u.DistributionMap(), 1, 1);

    compute_flux_theta_5(u, flux_t, geom);
    compute_flux_phi_5(u, flux_p, geom);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& ft = flux_t.const_array(mfi);
        auto const& fp = flux_p.const_array(mfi);
        auto const& rhs = Rhs.array(mfi);

        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double sin_t = std::sin(get_theta(i));

            double div_t = (ft(i,j,k,0) - ft(i-1,j,k,0)) / dtheta;
            double div_p = (fp(i,j,k,0) - fp(i,j-1,k,0)) / dphi;

            rhs(i,j,k) = -(div_t + div_p) / sin_t;
        });
    }
}

// ==============================================================================
// MAIN RUN LOOP
// ==============================================================================

int main(int argc, char* argv[]) {
    amrex::Initialize(argc, argv);
    {
        amrex::Print() << "Initializing AMReX Cosine Bell Advection (5th Order)...\n";
        amrex::Print() << "Resolution: " << N_THETA << " x " << N_PHI << "\n";

        amrex::Box domain(amrex::IntVect(AMREX_D_DECL(0,0,0)), amrex::IntVect(AMREX_D_DECL(N_THETA-1, N_PHI-1, 0)));

        amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(M_PI, 2.0*M_PI, 1.0)});
        amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 1, 0)};

        amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

        amrex::BoxArray ba(domain);
        ba.maxSize(amrex::IntVect(AMREX_D_DECL(16, N_PHI, 1)));
        amrex::DistributionMapping dm(ba);

        amrex::MultiFab u (ba, dm, 1, 3);
        amrex::MultiFab u0_mf(ba, dm, 1, 3);
        amrex::MultiFab u1(ba, dm, 1, 3);
        amrex::MultiFab u2(ba, dm, 1, 3);
        amrex::MultiFab u3(ba, dm, 1, 3);
        amrex::MultiFab u4(ba, dm, 1, 3);
        amrex::MultiFab rhs(ba, dm, 1, 0);
        amrex::MultiFab rhs_u3(ba, dm, 1, 0);

        for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
            auto const& arr = u.array(mfi);
            amrex::ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                arr(i,j,k) = exact_sol_bell(get_theta(i), get_phi(j));
            });
        }

        apply_boundaries(u, geom);

        double t = 0.0, t_end = 1.0;

        double theta_min = 0.5 * dtheta;
        double c_speed = u0 / dtheta + u0 / (std::sin(theta_min) * dphi);

        double base_CFL = 0.4;
        double dt = base_CFL / c_speed;

        amrex::Print() << "Starting Time Integration... CFL = " << base_CFL << " dt = " << dt << "\n";
        int step = 0;

        while (t < t_end) {
            if (t + dt > t_end) dt = t_end - t;

            amrex::MultiFab::Copy(u0_mf, u, 0, 0, 1, 3);

            // SSPRK Stage 1
            get_rhs_bell_5(u0_mf, rhs, geom);
            amrex::MultiFab::Copy(u1, u0_mf, 0, 0, 1, 0);
            amrex::MultiFab::Saxpy(u1, dt * 0.39175222700392, rhs, 0, 0, 1, 0);
            apply_boundaries(u1, geom);

            // Stage 2
            get_rhs_bell_5(u1, rhs, geom);
            amrex::MultiFab::LinComb(u2, 0.44437049406734, u0_mf, 0, 0.55562950593266, u1, 0, 0, 1, 0);
            amrex::MultiFab::Saxpy(u2, dt * 0.36841059262959, rhs, 0, 0, 1, 0);
            apply_boundaries(u2, geom);

            // Stage 3
            get_rhs_bell_5(u2, rhs, geom);
            amrex::MultiFab::LinComb(u3, 0.62010185138540, u0_mf, 0, 0.37989814861460, u2, 0, 0, 1, 0);
            amrex::MultiFab::Saxpy(u3, dt * 0.25189177424738, rhs, 0, 0, 1, 0);
            apply_boundaries(u3, geom);

            // Stage 4
            get_rhs_bell_5(u3, rhs_u3, geom);
            amrex::MultiFab::LinComb(u4, 0.17807995410773, u0_mf, 0, 0.82192004589227, u3, 0, 0, 1, 0);
            amrex::MultiFab::Saxpy(u4, dt * 0.54497475021237, rhs_u3, 0, 0, 1, 0);
            apply_boundaries(u4, geom);

            // Stage 5
            get_rhs_bell_5(u4, rhs, geom);
            u.setVal(0.0);
            amrex::MultiFab::Saxpy(u, 0.00683325884039, u0_mf, 0, 0, 1, 0);
            amrex::MultiFab::Saxpy(u, 0.51723167208978, u2, 0, 0, 1, 0);
            amrex::MultiFab::Saxpy(u, 0.12759831133288, u3, 0, 0, 1, 0);
            amrex::MultiFab::Saxpy(u, 0.34833675773694, u4, 0, 0, 1, 0);
            amrex::MultiFab::Saxpy(u, dt * 0.08460416338212, rhs_u3, 0, 0, 1, 0);
            amrex::MultiFab::Saxpy(u, dt * 0.22600748319395, rhs, 0, 0, 1, 0);
            apply_boundaries(u, geom);

            t += dt; step++;
            if (step % 100 == 0) amrex::Print() << "Step " << step << " | Time: " << t << " / " << t_end << "\n";
        }

        amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum,
                         amrex::ReduceOpSum, amrex::ReduceOpMax, amrex::ReduceOpMax> reduce_ops;
        amrex::ReduceData<double, double, double, double, double, double> reduce_data(reduce_ops);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
            auto const& arr = u.const_array(mfi);
            reduce_ops.eval(mfi.validbox(), reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
                double exact = exact_sol_bell(get_theta(i), get_phi(j));
                double err = std::abs(arr(i,j,k) - exact);

                double t_lo = i * dtheta;
                double t_hi = (i + 1) * dtheta;
                double vol = (std::cos(t_lo) - std::cos(t_hi)) * dphi;

                return {err * vol, std::abs(exact) * vol, err * err * vol, exact * exact * vol, err, std::abs(exact)};
            });
        }

        auto [num_l1, den_l1, num_l2, den_l2, num_li, den_li] = reduce_data.value();
        amrex::ParallelDescriptor::ReduceRealSum(num_l1);
        amrex::ParallelDescriptor::ReduceRealSum(den_l1);
        amrex::ParallelDescriptor::ReduceRealSum(num_l2);
        amrex::ParallelDescriptor::ReduceRealSum(den_l2);
        amrex::ParallelDescriptor::ReduceRealMax(num_li);
        amrex::ParallelDescriptor::ReduceRealMax(den_li);

        double l1_err = num_l1 / den_l1;
        double l2_err = std::sqrt(num_l2 / den_l2);
        double li_err = num_li / den_li;

        amrex::Print() << "\n--- Cosine Bell Advection Verification (5th Order) ---\n"
                       << std::scientific << std::setprecision(8)
                       << "l1 Error: " << l1_err << "\nl2 Error: " << l2_err
                       << "\nl-inf Error: " << li_err << "\nN_theta: " << N_THETA << " N_phi: " << N_PHI << "\n";

        amrex::WriteSingleLevelPlotfile("plt_final_bell", u, {"Q"}, geom, t, step);
    }
    amrex::Finalize();
    return 0;
}
