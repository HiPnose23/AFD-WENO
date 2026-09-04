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

// ==============================================================================
// CONFIG
// ==============================================================================
enum class Mode {RADIAL, MERIDIONAL, Z, AZIMUTHAL};


constexpr Mode TEST_MODE = Mode::MERIDIONAL;  
constexpr int ORDER = 5; // 3 or 5
constexpr int N = 64; 
constexpr double alpha_v = 1.0;

// RADIAL PARAMETERS
constexpr int m = 2;   // 0: Cartesian 1: Cylindrical 2: Spherical
constexpr double a_rad = 10.0;
constexpr double b_rad = 0.0; 

// MERIDIONAL PARAMETERS
// Only spherical
constexpr double a_mer = 10.0;
constexpr double b_mer = 0.0; 

// Z PARAMETERS
// Only cylindircal
constexpr double a_z = 5;
constexpr double b_z = 0.5;

// AZIMUTHAL PARAMETERS
constexpr double a_azi = 3.0;
constexpr double b_azi = M_PI;           
constexpr double r_fixed = 1.0;          // Fixed radius for the 1D ring
constexpr double theta_fixed = M_PI/2.0; // Fixed latitude (e.g., Equator)


// ==============================================================================
// RADIAL ADVECTION SOLVER
// ==============================================================================

const double dr = 2.0 / N; 
inline AMREX_GPU_DEVICE double get_r(int i) { return 0.0 + (i + 0.5) * dr; } 

inline AMREX_GPU_DEVICE double exact_sol_rad(double xi, double t) {
    double xi_0 = xi * std::exp(-alpha_v * t);
    double Q_0 = std::exp(-a_rad * a_rad * (xi_0 - b_rad) * (xi_0 - b_rad));
    return std::exp(-(m + 1.0) * alpha_v * t) * Q_0;
}

void apply_boundaries_rad(amrex::MultiFab& u, const amrex::Geometry& geom) {
    u.FillBoundary(geom.periodicity()); 

    const amrex::Box& domain = geom.Domain();
    int dom_lo_r = domain.smallEnd(0);
    int dom_hi_r = domain.bigEnd(0);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& arr = u.array(mfi);
        const amrex::Box& valid_box = mfi.validbox();
        amrex::Box grown_box = mfi.growntilebox();
        
        amrex::ParallelFor(grown_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            if (i < dom_lo_r) {
                int dist = dom_lo_r - i; 
                int sym_i = dom_lo_r + dist - 1; 
                arr(i, j, k) = arr(sym_i, j, k);
            } 
            else if (i > dom_hi_r) {
                arr(i, j, k) = arr(dom_hi_r, j, k);
            }
        });
    }
}

// ==============================================================================
// THIRD ORDER
// ==============================================================================

void compute_flux_rad(const amrex::MultiFab& u, amrex::MultiFab& flux, const amrex::Geometry& geom) {
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        amrex::Box flux_box = mfi.validbox(); 
        flux_box.growLo(0, 1); 
        
        auto const& u_arr = u.const_array(mfi);
        auto const& flx = flux.array(mfi);

        amrex::ParallelFor(flux_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double w_L, w_R, dwdr;
           
            auto point_f = [&](int idx) {
                double xi = get_r(idx);
                return std::pow(xi, m + 1) * alpha_v * u_arr(idx,j,k);
            };

            weno_ao_3_interpolation(
                u_arr(i-2,j,k), u_arr(i-1,j,k), u_arr(i,j,k), u_arr(i+1,j,k), u_arr(i+2,j,k),
                dr, w_L, w_R, dwdr);
            
            double f_face = w_R; 
            

            double d2f_r, d4f_r, d2dr2; // We only care about d2dr2 right now! 
            weno_ao_43_boundary(point_f(i-1), point_f(i), point_f(i+1), point_f(i+2), dr, d2f_r, d4f_r, d2dr2);

            double r_face = get_r(i) + 0.5 * dr;
            flx(i,j,k,0) = std::pow(r_face,m+1)*alpha_v*f_face - (dr*dr/24.0) * d2dr2;
        });
    }
}

void get_rhs_rad(const amrex::MultiFab& u, amrex::MultiFab& Rhs, const amrex::Geometry& geom) {
    amrex::MultiFab flux_act(u.boxArray(), u.DistributionMap(), 1, 1);
    compute_flux_rad(u, flux_act, geom);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto const& f_act = flux_act.const_array(mfi); 
        auto const& rhs   = Rhs.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double xi = get_r(i); 
            double div_F = (f_act(i,j,k,0) - f_act(i-1,j,k,0)) / (std::pow(xi, m) * dr);     
            rhs(i,j,k) = -div_F;
        });
    }
}

void run_radial_advection_3() {
    amrex::Print() << "Initializing AMReX Advection Solver (Radial)...\n";
    
    amrex::Box domain(amrex::IntVect(AMREX_D_DECL(0,0,0)), amrex::IntVect(AMREX_D_DECL(N-1,0,0)));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, -1.0, -1.0)}, {AMREX_D_DECL(2.0, 1.0, 1.0)});
    
    amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
    amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

    amrex::BoxArray ba(domain);
    ba.maxSize(16); 
    amrex::DistributionMapping dm(ba);

    amrex::MultiFab u(ba, dm, 1, 3);
    amrex::MultiFab u1(ba, dm, 1, 3);
    amrex::MultiFab u2(ba, dm, 1, 3);
    amrex::MultiFab rhs(ba, dm, 1, 0);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& arr = u.array(mfi);
        amrex::ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            arr(i,j,k) = exact_sol_rad(get_r(i), 0.0);
        });
    }
    apply_boundaries_rad(u, geom); 

    double t = 0.0, t_end = 1.0, CFL = 0.4;
    double c_speed = alpha_v * 2.0; 

    amrex::Print() << "Starting Time Integration...\n";
    BL_PROFILE_VAR("Evolution_Loop", pmain);
    int step = 0;
    double dt = CFL * dr / c_speed;

    while (t < t_end) {
        if (t + dt > t_end) dt = t_end - t;

        get_rhs_rad(u, rhs, geom);
        amrex::MultiFab::Copy(u1, u, 0, 0, 1, 0);                 
        amrex::MultiFab::Saxpy(u1, dt, rhs, 0, 0, 1, 0);          
        apply_boundaries_rad(u1, geom);

        get_rhs_rad(u1, rhs, geom);
        amrex::MultiFab::LinComb(u2, 0.75, u, 0, 0.25, u1, 0, 0, 1, 0); 
        amrex::MultiFab::Saxpy(u2, 0.25 * dt, rhs, 0, 0, 1, 0);
        apply_boundaries_rad(u2, geom);

        get_rhs_rad(u2, rhs, geom);
        amrex::MultiFab::LinComb(u, 1.0/3.0, u, 0, 2.0/3.0, u2, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u, (2.0/3.0) * dt, rhs, 0, 0, 1, 0);
        apply_boundaries_rad(u, geom);

        t += dt; step++;
        if (step % 10 == 0) amrex::Print() << "Step " << step << " | Time: " << t << " / " << t_end << "\n";
    }
    BL_PROFILE_VAR_STOP(pmain); 

    amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum> reduce_ops;
    amrex::ReduceData<double, double, double, double, double> reduce_data(reduce_ops);
    using ReduceTuple = typename amrex::ReduceData<double, double, double, double, double>::Type;

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
    auto const& arr = u.const_array(mfi);
    reduce_ops.eval(mfi.validbox(), reduce_data,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
            double exact = exact_sol_rad(get_r(i), t);
            double err   = std::abs(arr(i,j,k) - exact);

            double r_c = get_r(i);
            double w   = std::pow(r_c, (double)m) * dr;   // r^2 dr for m = 2

            return {err, err * w, err * err * w, err, w};
        });
    }

    ReduceTuple hv = reduce_data.value();
    double linf   = amrex::get<0>(hv);
    double num    = amrex::get<1>(hv);
    double l2num  = amrex::get<2>(hv);
    double sum    = amrex::get<3>(hv);
    double wsum   = amrex::get<4>(hv);

    amrex::ParallelDescriptor::ReduceRealMax(linf);
    amrex::ParallelDescriptor::ReduceRealSum(num);
    amrex::ParallelDescriptor::ReduceRealSum(l2num);
    amrex::ParallelDescriptor::ReduceRealSum(sum);
    amrex::ParallelDescriptor::ReduceRealSum(wsum);

    double l1_weighted = num / wsum;              // L1 variant for AFD
    double l1_plain    = sum / (double)N;    // unweighted
    double l2_weighted = std::sqrt(l2num / wsum);
    amrex::Print() << "\n--- Advection Benchmark Verification ---\n" << std::scientific << std::setprecision(8)
    << "L1 Error:    " << l1_plain << "\nL1 Error (Unweighted):    " << l1_plain << "\nL-inf Error: " << linf << "\nN: " << N << "\n";
    
    amrex::WriteSingleLevelPlotfile("plt_final", u, {"u"}, geom, t, step);
}

// ==============================================================================
// FIFTH ORDER
// ==============================================================================

void compute_flux_rad_5(const amrex::MultiFab& u, amrex::MultiFab& flux, const amrex::Geometry& geom) {
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        amrex::Box flux_box = mfi.validbox(); 
        flux_box.growLo(0, 1);  
        
        auto const& u_arr = u.const_array(mfi);
        auto const& flx = flux.array(mfi);

        amrex::ParallelFor(flux_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double w_L, w_R, dwdr;
           
            auto point_f = [&](int idx) {
                double xi = get_r(idx);
                return xi * alpha_v * u_arr(idx,j,k);
            };

            // 5th-order interpolation for the flux face value
            weno_ao_5_3_interpolation(
                u_arr(i-2,j,k), u_arr(i-1,j,k), u_arr(i,j,k), u_arr(i+1,j,k), u_arr(i+2,j,k),
                dr, w_L, w_R, dwdr);
            
            double f_face = w_R; 
            
            // 6th-order boundary interpolation for the exact 2nd and 4th derivatives
            double d2dr2, d4dr4; 
            weno_ao_63_boundary(
                point_f(i-2), point_f(i-1), point_f(i), 
                point_f(i+1), point_f(i+2), point_f(i+3), 
                dr, d2dr2, d4dr4);

            double r_face = get_r(i) + 0.5 * dr;
           
            
            flx(i,j,k,0) = r_face * alpha_v * f_face 
                         - (dr * dr / 24.0) * d2dr2 
                         + (7.0 * dr * dr * dr * dr / 5760.0) * d4dr4;
        });
    }
}

void get_rhs_rad_5(const amrex::MultiFab& u, amrex::MultiFab& Rhs, const amrex::Geometry& geom) {
    amrex::MultiFab flux_act(u.boxArray(), u.DistributionMap(), 1, 1);
    compute_flux_rad_5(u, flux_act, geom);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto const& f_act = flux_act.const_array(mfi); 
        auto const& rhs   = Rhs.array(mfi);
        auto const& u_arr = u.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double xi = get_r(i); 
            double div_F = (f_act(i,j,k,0) - f_act(i-1,j,k,0)) / dr;     
            double src = m * alpha_v * u_arr(i,j,k);
            rhs(i,j,k) = -div_F - src;
        });
    }
}

void run_radial_advection_5() {
    amrex::Print() << "Initializing AMReX Advection Solver (Radial 5th Order)...\n";
    
    amrex::Box domain(amrex::IntVect(AMREX_D_DECL(0,0,0)), amrex::IntVect(AMREX_D_DECL(N-1,0,0)));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, -1.0, -1.0)}, {AMREX_D_DECL(2.0, 1.0, 1.0)});
    
    amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
    amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

    amrex::BoxArray ba(domain);
    ba.maxSize(16); 
    amrex::DistributionMapping dm(ba);

    amrex::MultiFab u (ba, dm, 1, 3);
    amrex::MultiFab u0(ba, dm, 1, 3);
    amrex::MultiFab u1(ba, dm, 1, 3);
    amrex::MultiFab u2(ba, dm, 1, 3);
    amrex::MultiFab u3(ba, dm, 1, 3);
    amrex::MultiFab u4(ba, dm, 1, 3);
    amrex::MultiFab rhs(ba, dm, 1, 0);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& arr = u.array(mfi);
        amrex::ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            arr(i,j,k) = exact_sol_rad(get_r(i), 0.0);
        });
    }
    apply_boundaries_rad(u, geom); 

    double t = 0.0, t_end = 1.0;
    double c_speed = alpha_v * 2.0; 

    // Time scaling matching Section 6 of the Balsara paper 
    double N_base = 64.0;
    double base_CFL = 0.4;
    double scaled_CFL = base_CFL * std::pow(N_base / (double)N, 0.25); 
    double dt = scaled_CFL * dr / c_speed;

    amrex::Print() << "Starting Time Integration... Scaled CFL = " << scaled_CFL << "\n";
    BL_PROFILE_VAR("Evolution_Loop", pmain);
    int step = 0;

    while (t < t_end) {
        if (t + dt > t_end) dt = t_end - t;

        amrex::MultiFab::Copy(u0, u, 0, 0, 1, 3);

        // --- Stage 1 ---
        get_rhs_rad_5(u0, rhs, geom);
        amrex::MultiFab::Copy(u1, u0, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u1, dt * 0.39175222700392, rhs, 0, 0, 1, 0);
        apply_boundaries_rad(u1, geom);

        // --- Stage 2 ---
        get_rhs_rad_5(u1, rhs, geom);
        amrex::MultiFab::LinComb(u2, 0.44437049406734, u0, 0, 0.55562950593266, u1, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u2, dt * 0.36841059262959, rhs, 0, 0, 1, 0);
        apply_boundaries_rad(u2, geom);

        // --- Stage 3 ---
        get_rhs_rad_5(u2, rhs, geom);
        amrex::MultiFab::LinComb(u3, 0.62010185138540, u0, 0, 0.37989814861460, u2, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u3, dt * 0.25189177424738, rhs, 0, 0, 1, 0);
        apply_boundaries_rad(u3, geom);

        // --- Stage 4 ---
        amrex::MultiFab rhs_u3(ba, dm, 1, 0); // In order to save it and use it in stage 5!
        get_rhs_rad_5(u3, rhs_u3, geom);
        amrex::MultiFab::LinComb(u4, 0.17807995410773, u0, 0, 0.82192004589227, u3, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u4, dt * 0.54497475021237, rhs_u3, 0, 0, 1, 0);
        apply_boundaries_rad(u4, geom);

        // --- Stage 5  ---
        get_rhs_rad_5(u4, rhs, geom);
        u.setVal(0.0);
        amrex::MultiFab::Saxpy(u, 0.00683325884039, u0, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u, 0.51723167208978, u2, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u, 0.12759831133288, u3, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u, 0.34833675773694, u4, 0, 0, 1, 0);

        amrex::MultiFab::Saxpy(u, dt * 0.08460416338212, rhs_u3, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u, dt * 0.22600748319395, rhs, 0, 0, 1, 0);   
        apply_boundaries_rad(u, geom);

        t += dt; 
        step++;
        if (step % 10 == 0) amrex::Print() << "Step " << step << " | Time: " << t << " / " << t_end << "\n";
    }
    BL_PROFILE_VAR_STOP(pmain); 

    amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum> reduce_ops;
    amrex::ReduceData<double, double, double, double, double> reduce_data(reduce_ops);
    using ReduceTuple = typename amrex::ReduceData<double, double, double, double, double>::Type;

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
    auto const& arr = u.const_array(mfi);
    reduce_ops.eval(mfi.validbox(), reduce_data,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
            double exact = exact_sol_rad(get_r(i), t);
            double err   = std::abs(arr(i,j,k) - exact);

            double r_c = get_r(i);
            double w   = std::pow(r_c, (double)m) * dr;   // r^2 dr for m = 2

            return {err, err * w, err * err * w, err, w};
        });
    }

    ReduceTuple hv = reduce_data.value();
    double linf   = amrex::get<0>(hv);
    double num    = amrex::get<1>(hv);
    double l2num  = amrex::get<2>(hv);
    double sum    = amrex::get<3>(hv);
    double wsum   = amrex::get<4>(hv);

    amrex::ParallelDescriptor::ReduceRealMax(linf);
    amrex::ParallelDescriptor::ReduceRealSum(num);
    amrex::ParallelDescriptor::ReduceRealSum(l2num);
    amrex::ParallelDescriptor::ReduceRealSum(sum);
    amrex::ParallelDescriptor::ReduceRealSum(wsum);

    double l1_weighted = num / wsum;              // L1 variant for AFD
    double l1_plain    = sum / (double)N;    // unweighted
    double l2_weighted = std::sqrt(l2num / wsum);
    amrex::Print() << "\n--- Advection Benchmark Verification ---\n" << std::scientific << std::setprecision(8)
    << "L1 Error:    " << l1_plain << "\nL1 Error (Unweighted):    " << l1_plain << "\nL-inf Error: " << linf << "\nN: " << N << "\n";
    
    amrex::WriteSingleLevelPlotfile("plt_final", u, {"u"}, geom, t, step);
}

// ==============================================================================
// MERIDIONAL ADVECTION SOLVER
// ==============================================================================

// ==============================================================================
// THIRD ORDER
// ==============================================================================

const double dtheta = (M_PI / 2.0) / N; 
inline AMREX_GPU_DEVICE double get_theta(int i) { return 0.0 + (i + 0.5) * dtheta; } 

inline AMREX_GPU_DEVICE double exact_sol_mer(double theta, double t) {
    double theta_0 = theta * std::exp(-alpha_v * t);
    double Q_0 = 0.0;
    if (std::abs(theta_0 - b_mer) < (M_PI / a_mer)) {
        double arg = (1.0 + std::cos(a_mer * (theta_0 - b_mer))) / 2.0;
        Q_0 = arg * arg * arg * arg;
    }
    double exp = std::sin(theta_0) / std::sin(theta);
    return std::exp(-alpha_v * t) * exp * Q_0;
}

void apply_boundaries_mer(amrex::MultiFab& u, const amrex::Geometry& geom) {
    u.FillBoundary(geom.periodicity()); 

    const amrex::Box& domain = geom.Domain();
    int dom_lo_r = domain.smallEnd(0);
    int dom_hi_r = domain.bigEnd(0);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& arr = u.array(mfi);
        const amrex::Box& valid_box = mfi.validbox();
        amrex::Box grown_box = mfi.growntilebox();
        
        amrex::ParallelFor(grown_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            if (i < dom_lo_r) {
                int dist = dom_lo_r - i; 
                int sym_i = dom_lo_r + dist - 1; 
                arr(i, j, k) = arr(sym_i, j, k);
            } 
            else if (i > dom_hi_r) {
                arr(i, j, k) = arr(dom_hi_r, j, k);
            }
        });
    }
}

void compute_flux_mer(const amrex::MultiFab& u, amrex::MultiFab& flux, const amrex::Geometry& geom) {
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        amrex::Box flux_box = mfi.validbox(); 
        flux_box.growLo(0, 1); 
        
        auto const& u_arr = u.const_array(mfi);
        auto const& flx = flux.array(mfi);

        amrex::ParallelFor(flux_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double w_L, w_R, dwdr;
           
            auto point_f = [&](int idx) {
                double theta = get_theta(idx);
                return std::sin(theta) * (alpha_v * theta) * u_arr(idx,j,k);            
            };

            weno_ao_3_interpolation(
                u_arr(i-2,j,k), u_arr(i-1,j,k), u_arr(i,j,k), u_arr(i+1,j,k), u_arr(i+2,j,k),
                dtheta, w_L, w_R, dwdr);
            
            double f_face = w_R; 
            
            double d2f_theta, d4f_theta, d2dtheta2; // We only care about d2dtheta2 right now! 
            weno_ao_43_boundary(point_f(i-1), point_f(i), point_f(i+1), point_f(i+2), dtheta, d2f_theta, d4f_theta, d2dtheta2);

            
            double theta_face = get_theta(i) + 0.5 * dtheta;
            flx(i,j,k,0) = std::sin(theta_face) * (alpha_v * theta_face) * f_face - (dtheta*dtheta/24.0) * d2dtheta2;
        });
    }
}

void get_rhs_mer(const amrex::MultiFab& u, amrex::MultiFab& Rhs, const amrex::Geometry& geom) {
    amrex::MultiFab flux_act(u.boxArray(), u.DistributionMap(), 1, 1);
    compute_flux_mer(u, flux_act, geom);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto const& f_act = flux_act.const_array(mfi); 
        auto const& rhs   = Rhs.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double theta = get_theta(i); 
            double div_F = (f_act(i,j,k,0) - f_act(i-1,j,k,0)) / dtheta;     
            rhs(i,j,k) = -(1.0 / std::sin(theta)) * div_F;
        });
    }
}

void run_meridional_advection_3() {
    amrex::Print() << "Initializing AMReX Advection Solver (Meridional)...\n";
    
    amrex::Box domain(amrex::IntVect(AMREX_D_DECL(0,0,0)), amrex::IntVect(AMREX_D_DECL(N-1,0,0)));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, -1.0, -1.0)}, {AMREX_D_DECL(M_PI/2.0, 1.0, 1.0)});
    
    amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
    amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

    amrex::BoxArray ba(domain);
    ba.maxSize(16); 
    amrex::DistributionMapping dm(ba);

    amrex::MultiFab u(ba, dm, 1, 3);
    amrex::MultiFab u1(ba, dm, 1, 3);
    amrex::MultiFab u2(ba, dm, 1, 3);
    amrex::MultiFab rhs(ba, dm, 1, 0);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& arr = u.array(mfi);
        amrex::ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            arr(i,j,k) = exact_sol_mer(get_theta(i), 0.0);
        });
    }
    apply_boundaries_mer(u, geom); 

    double t = 0.0, t_end = 1.0, CFL = 0.4;
    double c_speed = alpha_v * (M_PI / 2.0); 

    amrex::Print() << "Starting Time Integration...\n";
    BL_PROFILE_VAR("Evolution_Loop", pmain);
    int step = 0;
    double dt = CFL * dtheta / c_speed;

    while (t < t_end) {
        if (t + dt > t_end) dt = t_end - t;

        get_rhs_mer(u, rhs, geom);
        amrex::MultiFab::Copy(u1, u, 0, 0, 1, 0);                 
        amrex::MultiFab::Saxpy(u1, dt, rhs, 0, 0, 1, 0);          
        apply_boundaries_mer(u1, geom);

        get_rhs_mer(u1, rhs, geom);
        amrex::MultiFab::LinComb(u2, 0.75, u, 0, 0.25, u1, 0, 0, 1, 0); 
        amrex::MultiFab::Saxpy(u2, 0.25 * dt, rhs, 0, 0, 1, 0);
        apply_boundaries_mer(u2, geom);

        get_rhs_mer(u2, rhs, geom);
        amrex::MultiFab::LinComb(u, 1.0/3.0, u, 0, 2.0/3.0, u2, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u, (2.0/3.0) * dt, rhs, 0, 0, 1, 0);
        apply_boundaries_mer(u, geom);

        t += dt; step++;
        if (step % 10 == 0) amrex::Print() << "Step " << step << " | Time: " << t << " / " << t_end << "\n";
    }
    BL_PROFILE_VAR_STOP(pmain); 

   
 amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum> reduce_ops;
    amrex::ReduceData<double, double, double, double, double> reduce_data(reduce_ops);
    using ReduceTuple = typename amrex::ReduceData<double, double, double, double, double>::Type;

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
    auto const& arr = u.const_array(mfi);
    reduce_ops.eval(mfi.validbox(), reduce_data,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
            double exact = exact_sol_mer(get_theta(i), t);
            double err   = std::abs(arr(i,j,k) - exact);

            double r_c = get_r(i);
            double w   = std::pow(r_c, (double)m) * dr;   // r^2 dr for m = 2

            return {err, err * w, err * err * w, err, w};
        });
    }
    // Try new L1 which is like the inf but reducesum + divide by N

    ReduceTuple hv = reduce_data.value();
    double linf   = amrex::get<0>(hv);
    double num    = amrex::get<1>(hv);
    double l2num  = amrex::get<2>(hv);
    double sum    = amrex::get<3>(hv);
    double wsum   = amrex::get<4>(hv);

    amrex::ParallelDescriptor::ReduceRealMax(linf);
    amrex::ParallelDescriptor::ReduceRealSum(num);
    amrex::ParallelDescriptor::ReduceRealSum(l2num);
    amrex::ParallelDescriptor::ReduceRealSum(sum);
    amrex::ParallelDescriptor::ReduceRealSum(wsum);

    double l1_weighted = num / wsum;              // L1 variant for AFD
    double l1_plain    = sum / (double)N;    // unweighted
    double l2_weighted = std::sqrt(l2num / wsum);
    amrex::Print() << "\n--- Advection Benchmark Verification ---\n" << std::scientific << std::setprecision(8)
    << "L1 Error:    " << l1_plain << "\nL1 Error (Unweighted):    " << l1_plain << "\nL-inf Error: " << linf << "\nN: " << N << "\n";
    
    amrex::WriteSingleLevelPlotfile("plt_final", u, {"u"}, geom, t, step);
}

// ==============================================================================
// FIFTH ORDER
// ==============================================================================

void compute_flux_mer_5(const amrex::MultiFab& u, amrex::MultiFab& flux, const amrex::Geometry& geom) {
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        amrex::Box flux_box = mfi.validbox(); 
        flux_box.growLo(0, 1); 
        
        auto const& u_arr = u.const_array(mfi);
        auto const& flx = flux.array(mfi);

        amrex::ParallelFor(flux_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double w_L, w_R, dwdr;
           
            auto point_f = [&](int idx) {
                double theta = get_theta(idx);
                return std::sin(theta) * (alpha_v * theta) * u_arr(idx,j,k);            
            };

            weno_ao_5_3_interpolation(
                u_arr(i-2,j,k), u_arr(i-1,j,k), u_arr(i,j,k), u_arr(i+1,j,k), u_arr(i+2,j,k),
                dtheta, w_L, w_R, dwdr);
            
            double f_face = w_R; 
            
            double d2dtheta2, d4dtheta4; 
            weno_ao_63_boundary(
                point_f(i-2), point_f(i-1), point_f(i), 
                point_f(i+1), point_f(i+2), point_f(i+3), 
                dtheta, d2dtheta2, d4dtheta4);

            double theta_face = get_theta(i) + 0.5 * dtheta;
            
            flx(i,j,k,0) = std::sin(theta_face) * (alpha_v * theta_face) * f_face 
                         - (dtheta * dtheta / 24.0) * d2dtheta2
                         + (7.0 * dtheta * dtheta * dtheta * dtheta / 5760.0) * d4dtheta4;
        });
    }
}

void get_rhs_mer_5(const amrex::MultiFab& u, amrex::MultiFab& Rhs, const amrex::Geometry& geom) {
    amrex::MultiFab flux_act(u.boxArray(), u.DistributionMap(), 1, 1);
    compute_flux_mer_5(u, flux_act, geom);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto const& f_act = flux_act.const_array(mfi); 
        auto const& rhs   = Rhs.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double theta = get_theta(i); 
            double div_F = (f_act(i,j,k,0) - f_act(i-1,j,k,0)) / dtheta;     
            rhs(i,j,k) = -(1.0 / std::sin(theta)) * div_F;
        });
    }
}

void run_meridional_advection_5() {
    amrex::Print() << "Initializing AMReX Advection Solver (Meridional 5th Order)...\n";
    
    amrex::Box domain(amrex::IntVect(AMREX_D_DECL(0,0,0)), amrex::IntVect(AMREX_D_DECL(N-1,0,0)));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, -1.0, -1.0)}, {AMREX_D_DECL(M_PI/2.0, 1.0, 1.0)});
    
    amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
    amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

    amrex::BoxArray ba(domain);
    ba.maxSize(16); 
    amrex::DistributionMapping dm(ba);

    amrex::MultiFab u (ba, dm, 1, 3);
    amrex::MultiFab u0(ba, dm, 1, 3);
    amrex::MultiFab u1(ba, dm, 1, 3);
    amrex::MultiFab u2(ba, dm, 1, 3);
    amrex::MultiFab u3(ba, dm, 1, 3);
    amrex::MultiFab u4(ba, dm, 1, 3);
    amrex::MultiFab rhs(ba, dm, 1, 0);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& arr = u.array(mfi);
        amrex::ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            arr(i,j,k) = exact_sol_mer(get_theta(i), 0.0);
        });
    }
    apply_boundaries_mer(u, geom); 

    double t = 0.0, t_end = 1.0;
    double c_speed = alpha_v * (M_PI / 2.0); 

    double N_base = 64.0;
    double base_CFL = 0.4;
    double scaled_CFL = base_CFL * std::pow(N_base / (double)N, 0.25);
    double dt = scaled_CFL * dtheta / c_speed;

    amrex::Print() << "Starting Time Integration... Scaled CFL = " << scaled_CFL << "\n";
    BL_PROFILE_VAR("Evolution_Loop", pmain);
    int step = 0;

    while (t < t_end) {
        if (t + dt > t_end) dt = t_end - t;

        amrex::MultiFab::Copy(u0, u, 0, 0, 1, 3);

        // --- Stage 1 ---
        get_rhs_mer_5(u0, rhs, geom);
        amrex::MultiFab::Copy(u1, u0, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u1, dt * 0.39175222700392, rhs, 0, 0, 1, 0);
        apply_boundaries_mer(u1, geom);

        // --- Stage 2 ---
        get_rhs_mer_5(u1, rhs, geom);
        amrex::MultiFab::LinComb(u2, 0.44437049406734, u0, 0, 0.55562950593266, u1, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u2, dt * 0.36841059262959, rhs, 0, 0, 1, 0);
        apply_boundaries_mer(u2, geom);

        // --- Stage 3 ---
        get_rhs_mer_5(u2, rhs, geom);
        amrex::MultiFab::LinComb(u3, 0.62010185138540, u0, 0, 0.37989814861460, u2, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u3, dt * 0.25189177424738, rhs, 0, 0, 1, 0);
        apply_boundaries_mer(u3, geom);

        // --- Stage 4 ---
        amrex::MultiFab rhs_u3(ba, dm, 1, 0); 
        get_rhs_mer_5(u3, rhs_u3, geom);
        amrex::MultiFab::LinComb(u4, 0.17807995410773, u0, 0, 0.82192004589227, u3, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u4, dt * 0.54497475021237, rhs_u3, 0, 0, 1, 0);
        apply_boundaries_mer(u4, geom);

        // --- Stage 5 ---
        get_rhs_mer_5(u4, rhs, geom);
        u.setVal(0.0);
        amrex::MultiFab::Saxpy(u, 0.00683325884039, u0, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u, 0.51723167208978, u2, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u, 0.12759831133288, u3, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u, 0.34833675773694, u4, 0, 0, 1, 0);

        amrex::MultiFab::Saxpy(u, dt * 0.08460416338212, rhs_u3, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u, dt * 0.22600748319395, rhs, 0, 0, 1, 0);   
        apply_boundaries_mer(u, geom);

        t += dt; 
        step++;
        if (step % 10 == 0) amrex::Print() << "Step " << step << " | Time: " << t << " / " << t_end << "\n";
    }
    BL_PROFILE_VAR_STOP(pmain); 

    
 amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum> reduce_ops;
    amrex::ReduceData<double, double, double, double, double> reduce_data(reduce_ops);
    using ReduceTuple = typename amrex::ReduceData<double, double, double, double, double>::Type;

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
    auto const& arr = u.const_array(mfi);
    reduce_ops.eval(mfi.validbox(), reduce_data,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
            double exact = exact_sol_mer(get_theta(i), t);
            double err   = std::abs(arr(i,j,k) - exact);

            double r_c = get_r(i);
            double w   = std::pow(r_c, (double)m) * dr;   // r^2 dr for m = 2

            return {err, err * w, err * err * w, err, w};
        });
    }


    ReduceTuple hv = reduce_data.value();
    double linf   = amrex::get<0>(hv);
    double num    = amrex::get<1>(hv);
    double l2num  = amrex::get<2>(hv);
    double sum    = amrex::get<3>(hv);
    double wsum   = amrex::get<4>(hv);

    amrex::ParallelDescriptor::ReduceRealMax(linf);
    amrex::ParallelDescriptor::ReduceRealSum(num);
    amrex::ParallelDescriptor::ReduceRealSum(l2num);
    amrex::ParallelDescriptor::ReduceRealSum(sum);
    amrex::ParallelDescriptor::ReduceRealSum(wsum);

    double l1_weighted = num / wsum;              // L1 variant for AFD
    double l1_plain    = sum / (double)N;    // unweighted
    double l2_weighted = std::sqrt(l2num / wsum);
    amrex::Print() << "\n--- Advection Benchmark Verification ---\n" << std::scientific << std::setprecision(8)
    << "L1 Error:    " << l1_plain << "\nL1 Error (Unweighted):    " << l1_plain << "\nL-inf Error: " << linf << "\nN: " << N << "\n";
    
    amrex::WriteSingleLevelPlotfile("plt_final", u, {"u"}, geom, t, step);
}

// ==============================================================================
// 3. Z-CARTESIAN ADVECTION SOLVER
// ==============================================================================

// ==============================================================================
// THIRD ORDER
// ==============================================================================

const double dz = 2.0 / N; 
inline AMREX_GPU_DEVICE double get_z(int i) { return -1.0 + (i + 0.5) * dz; } 

inline AMREX_GPU_DEVICE double exact_sol_z(double z, double t) {
    double z_0 = z - alpha_v * t;
    // Periodic wrap around if it leaves [-1, 1]
    while (z_0 > 1.0) z_0 -= 2.0;
    while (z_0 < -1.0) z_0 += 2.0;
    return std::exp(-a_z * a_z * (z_0 - b_z) * (z_0 - b_z));
}

void apply_boundaries_z(amrex::MultiFab& u, const amrex::Geometry& geom) {
    u.FillBoundary(geom.periodicity()); 
}

void compute_flux_z(const amrex::MultiFab& u, amrex::MultiFab& flux, const amrex::Geometry& geom) {
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        amrex::Box flux_box = mfi.validbox(); flux_box.growLo(0, 1); 

        auto const& u_arr = u.const_array(mfi); auto const& flx = flux.array(mfi);

        amrex::ParallelFor(flux_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double w_L, w_R, dwdz;

            auto point_f = [&](int idx) { 
                return alpha_v * u_arr(idx,j,k);
            };
            
            weno_ao_3_interpolation(
                u_arr(i-2,j,k), u_arr(i-1,j,k), u_arr(i,j,k), u_arr(i+1,j,k), u_arr(i+2,j,k),
                dz, w_L, w_R, dwdz);
            
            double f_face = w_R;
            double d2f_z, d4f_z, d2dz2; // We only care about d2dz2 right now! 
            weno_ao_43_boundary(point_f(i-1), point_f(i), point_f(i+1), point_f(i+2), dz, d2f_z, d4f_z, d2dz2);

                            
            flx(i,j,k,0) = alpha_v * f_face - (dz*dz/24.0) * d2dz2;
        });
    }
}

void get_rhs_z(const amrex::MultiFab& u, amrex::MultiFab& Rhs, const amrex::Geometry& geom) {
    amrex::MultiFab flux_act(u.boxArray(), u.DistributionMap(), 1, 1); 
    compute_flux_z(u, flux_act, geom);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& f_act = flux_act.const_array(mfi);
        auto const& rhs = Rhs.array(mfi);
        
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            rhs(i,j,k) = -((f_act(i,j,k,0) - f_act(i-1,j,k,0)) / dz); 
        });
    }
}

void run_z_advection_3() {
    amrex::Print() << "Initializing AMReX Advection Solver (Z-Cartesian)...\n";
    
    amrex::Box domain(amrex::IntVect(AMREX_D_DECL(0,0,0)), amrex::IntVect(AMREX_D_DECL(N-1,0,0)));
    amrex::RealBox real_box({AMREX_D_DECL(-1.0, -1.0, -1.0)}, {AMREX_D_DECL(1.0, 1.0, 1.0)});
    
    amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 0, 0)};
    amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

    amrex::BoxArray ba(domain);
    ba.maxSize(16);
    amrex::DistributionMapping dm(ba);

    amrex::MultiFab u(ba, dm, 1, 3);
    amrex::MultiFab u1(ba, dm, 1, 3);
    amrex::MultiFab u2(ba, dm, 1, 3);
    amrex::MultiFab rhs(ba, dm, 1, 0);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& arr = u.array(mfi);
        amrex::ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                arr(i,j,k) = exact_sol_z(get_z(i), 0.0);
        });
    }
    apply_boundaries_z(u, geom); 

    double t = 0.0, t_end = 1.0, CFL = 0.4;
    double c_speed = alpha_v;

    amrex::Print() << "Starting Time Integration...\n";
    BL_PROFILE_VAR("Evolution_Loop", pmain);
    int step = 0;
    double dt = CFL * dz / c_speed;

    while (t < t_end) {
        if (t + dt > t_end) dt = t_end - t;

        get_rhs_z(u, rhs, geom);
        amrex::MultiFab::Copy(u1, u, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u1, dt, rhs, 0, 0, 1, 0);
        apply_boundaries_z(u1, geom);

        get_rhs_z(u1, rhs, geom);
        amrex::MultiFab::LinComb(u2, 0.75, u, 0, 0.25, u1, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u2, 0.25 * dt, rhs, 0, 0, 1, 0);
        apply_boundaries_z(u2, geom);

        get_rhs_z(u2, rhs, geom);
        amrex::MultiFab::LinComb(u, 1.0/3.0, u, 0, 2.0/3.0, u2, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u, (2.0/3.0) * dt, rhs, 0, 0, 1, 0);
        apply_boundaries_z(u, geom);
        
        t += dt; step++;
        if (step % 10 == 0) amrex::Print() << "Step " << step << " | Time: " << t << " / " << t_end << "\n";

    }
    BL_PROFILE_VAR_STOP(pmain); 

    
    amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum> reduce_ops;
    amrex::ReduceData<double, double, double, double, double> reduce_data(reduce_ops);
    using ReduceTuple = typename amrex::ReduceData<double, double, double, double, double>::Type;

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
    auto const& arr = u.const_array(mfi);
    reduce_ops.eval(mfi.validbox(), reduce_data,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
            double exact = exact_sol_z(get_z(i), t);
            double err   = std::abs(arr(i,j,k) - exact);

            double r_c = get_r(i);
            double w   = std::pow(r_c, (double)m) * dr;   // r^2 dr for m = 2

            return {err, err * w, err * err * w, err, w};
        });
    }

    ReduceTuple hv = reduce_data.value();
    double linf   = amrex::get<0>(hv);
    double num    = amrex::get<1>(hv);
    double l2num  = amrex::get<2>(hv);
    double sum    = amrex::get<3>(hv);
    double wsum   = amrex::get<4>(hv);

    amrex::ParallelDescriptor::ReduceRealMax(linf);
    amrex::ParallelDescriptor::ReduceRealSum(num);
    amrex::ParallelDescriptor::ReduceRealSum(l2num);
    amrex::ParallelDescriptor::ReduceRealSum(sum);
    amrex::ParallelDescriptor::ReduceRealSum(wsum);

    double l1_weighted = num / wsum;              // L1 variant for AFD
    double l1_plain    = sum / (double)N;    // unweighted
    double l2_weighted = std::sqrt(l2num / wsum);
    amrex::Print() << "\n--- Advection Benchmark Verification ---\n" << std::scientific << std::setprecision(8)
    << "L1 Error:    " << l1_plain << "\nL1 Error (Unweighted):    " << l1_plain << "\nL-inf Error: " << linf << "\nN: " << N << "\n";
    
    amrex::WriteSingleLevelPlotfile("plt_final", u, {"u"}, geom, t, step);
}

// ==============================================================================
// FIFTH ORDER
// ==============================================================================

void compute_flux_z_5(const amrex::MultiFab& u, amrex::MultiFab& flux, const amrex::Geometry& geom) {
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        amrex::Box flux_box = mfi.validbox(); flux_box.growLo(0, 1); 

        auto const& u_arr = u.const_array(mfi); auto const& flx = flux.array(mfi);

        amrex::ParallelFor(flux_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double w_L, w_R, dwdz;

            auto point_f = [&](int idx) { 
                return alpha_v * u_arr(idx,j,k);
            };
            
            weno_ao_5_3_interpolation(
                u_arr(i-2,j,k), u_arr(i-1,j,k), u_arr(i,j,k), u_arr(i+1,j,k), u_arr(i+2,j,k),
                dz, w_L, w_R, dwdz);
            
            double f_face = w_R;
            
            double d2dz2, d4dz4; 
            weno_ao_63_boundary(
                point_f(i-2), point_f(i-1), point_f(i), 
                point_f(i+1), point_f(i+2), point_f(i+3), 
                dz, d2dz2, d4dz4);

            flx(i,j,k,0) = alpha_v * f_face 
                         - (dz * dz / 24.0) * d2dz2
                         + (7.0 * dz * dz * dz * dz / 5760.0) * d4dz4;
        });
    }
}

void get_rhs_z_5(const amrex::MultiFab& u, amrex::MultiFab& Rhs, const amrex::Geometry& geom) {
    amrex::MultiFab flux_act(u.boxArray(), u.DistributionMap(), 1, 1); 
    compute_flux_z_5(u, flux_act, geom);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& f_act = flux_act.const_array(mfi);
        auto const& rhs = Rhs.array(mfi);
        
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            rhs(i,j,k) = -((f_act(i,j,k,0) - f_act(i-1,j,k,0)) / dz); 
        });
    }
}

void run_z_advection_5() {
    amrex::Print() << "Initializing AMReX Advection Solver (Z-Cartesian 5th Order)...\n";
    
    amrex::Box domain(amrex::IntVect(AMREX_D_DECL(0,0,0)), amrex::IntVect(AMREX_D_DECL(N-1,0,0)));
    amrex::RealBox real_box({AMREX_D_DECL(-1.0, -1.0, -1.0)}, {AMREX_D_DECL(1.0, 1.0, 1.0)});
    
    amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 0, 0)};
    amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

    amrex::BoxArray ba(domain);
    ba.maxSize(16);
    amrex::DistributionMapping dm(ba);

    amrex::MultiFab u (ba, dm, 1, 3);
    amrex::MultiFab u0(ba, dm, 1, 3);
    amrex::MultiFab u1(ba, dm, 1, 3);
    amrex::MultiFab u2(ba, dm, 1, 3);
    amrex::MultiFab u3(ba, dm, 1, 3);
    amrex::MultiFab u4(ba, dm, 1, 3);
    amrex::MultiFab rhs(ba, dm, 1, 0);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& arr = u.array(mfi);
        amrex::ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                arr(i,j,k) = exact_sol_z(get_z(i), 0.0);
        });
    }
    apply_boundaries_z(u, geom); 

    double t = 0.0, t_end = 1.0;
    double c_speed = alpha_v;

    double N_base = 64.0;
    double base_CFL = 0.4;
    double scaled_CFL = base_CFL * std::pow(N_base / (double)N, 0.25);
    double dt = scaled_CFL * dz / c_speed;

    amrex::Print() << "Starting Time Integration... Scaled CFL = " << scaled_CFL << "\n";
    BL_PROFILE_VAR("Evolution_Loop", pmain);
    int step = 0;

    while (t < t_end) {
        if (t + dt > t_end) dt = t_end - t;

        amrex::MultiFab::Copy(u0, u, 0, 0, 1, 3);

        // --- Stage 1 ---
        get_rhs_z_5(u0, rhs, geom);
        amrex::MultiFab::Copy(u1, u0, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u1, dt * 0.39175222700392, rhs, 0, 0, 1, 0);
        apply_boundaries_z(u1, geom);

        // --- Stage 2 ---
        get_rhs_z_5(u1, rhs, geom);
        amrex::MultiFab::LinComb(u2, 0.44437049406734, u0, 0, 0.55562950593266, u1, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u2, dt * 0.36841059262959, rhs, 0, 0, 1, 0);
        apply_boundaries_z(u2, geom);

        // --- Stage 3 ---
        get_rhs_z_5(u2, rhs, geom);
        amrex::MultiFab::LinComb(u3, 0.62010185138540, u0, 0, 0.37989814861460, u2, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u3, dt * 0.25189177424738, rhs, 0, 0, 1, 0);
        apply_boundaries_z(u3, geom);

        // --- Stage 4 ---
        amrex::MultiFab rhs_u3(ba, dm, 1, 0); 
        get_rhs_z_5(u3, rhs_u3, geom);
        amrex::MultiFab::LinComb(u4, 0.17807995410773, u0, 0, 0.82192004589227, u3, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u4, dt * 0.54497475021237, rhs_u3, 0, 0, 1, 0);
        apply_boundaries_z(u4, geom);

        // --- Stage 5 ---
        get_rhs_z_5(u4, rhs, geom);
        u.setVal(0.0);
        amrex::MultiFab::Saxpy(u, 0.00683325884039, u0, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u, 0.51723167208978, u2, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u, 0.12759831133288, u3, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u, 0.34833675773694, u4, 0, 0, 1, 0);

        amrex::MultiFab::Saxpy(u, dt * 0.08460416338212, rhs_u3, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u, dt * 0.22600748319395, rhs, 0, 0, 1, 0);   
        apply_boundaries_z(u, geom);

        t += dt; 
        step++;
        if (step % 10 == 0) amrex::Print() << "Step " << step << " | Time: " << t << " / " << t_end << "\n";
    }
    BL_PROFILE_VAR_STOP(pmain); 

    
    amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum> reduce_ops;
    amrex::ReduceData<double, double, double, double, double> reduce_data(reduce_ops);
    using ReduceTuple = typename amrex::ReduceData<double, double, double, double, double>::Type;

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
    auto const& arr = u.const_array(mfi);
    reduce_ops.eval(mfi.validbox(), reduce_data,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
            double exact = exact_sol_z(get_z(i), t);
            double err   = std::abs(arr(i,j,k) - exact);

            double r_c = get_r(i);
            double w   = std::pow(r_c, (double)m) * dr;   // r^2 dr for m = 2

            return {err, err * w, err * err * w, err, w};
        });
    }

    ReduceTuple hv = reduce_data.value();
    double linf   = amrex::get<0>(hv);
    double num    = amrex::get<1>(hv);
    double l2num  = amrex::get<2>(hv);
    double sum    = amrex::get<3>(hv);
    double wsum   = amrex::get<4>(hv);

    amrex::ParallelDescriptor::ReduceRealMax(linf);
    amrex::ParallelDescriptor::ReduceRealSum(num);
    amrex::ParallelDescriptor::ReduceRealSum(l2num);
    amrex::ParallelDescriptor::ReduceRealSum(sum);
    amrex::ParallelDescriptor::ReduceRealSum(wsum);

    double l1_weighted = num / wsum;              // L1 variant for AFD
    double l1_plain    = sum / (double)N;    // unweighted
    double l2_weighted = std::sqrt(l2num / wsum);
    amrex::Print() << "\n--- Advection Benchmark Verification ---\n" << std::scientific << std::setprecision(8)
    << "L1 Error:    " << l1_plain << "\nL1 Error (Unweighted):    " << l1_plain << "\nL-inf Error: " << linf << "\nN: " << N << "\n";
    
    amrex::WriteSingleLevelPlotfile("plt_final", u, {"u"}, geom, t, step);
}

// ==============================================================================
// AZIMUTHAL ADVECTION SOLVER
// ==============================================================================

// ==============================================================================
// THIRD ORDER
// ==============================================================================

const double dphi = (2.0 * M_PI) / N; 
inline AMREX_GPU_DEVICE double get_phi(int i) { return 0.0 + (i + 0.5) * dphi; } 

inline AMREX_GPU_DEVICE double exact_sol_azi(double phi, double t) {
    double factor = (m == 1) ? (1.0 / r_fixed) : (1.0 / (r_fixed * std::sin(theta_fixed)));
    double omega = alpha_v * factor;
    double phi_0 = phi - omega * t;
    
    // Periodic wrap around for [0, 2*PI]
    while (phi_0 > 2.0 * M_PI) phi_0 -= 2.0 * M_PI;
    while (phi_0 < 0.0) phi_0 += 2.0 * M_PI;
    
    return std::exp(-a_azi * a_azi * (phi_0 - b_azi) * (phi_0 - b_azi));
}

void apply_boundaries_azi(amrex::MultiFab& u, const amrex::Geometry& geom) {
    u.FillBoundary(geom.periodicity()); 
}

void compute_flux_azi(const amrex::MultiFab& u, amrex::MultiFab& flux, const amrex::Geometry& geom) {
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        amrex::Box flux_box = mfi.validbox(); flux_box.growLo(0, 1); 

        auto const& u_arr = u.const_array(mfi); 
        auto const& flx = flux.array(mfi);

        amrex::ParallelFor(flux_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double w_L, w_R, dwdphi;

            auto point_f = [&](int idx) { 
                return alpha_v * u_arr(idx,j,k); 
            };
            
            weno_ao_3_interpolation(
                u_arr(i-2,j,k), u_arr(i-1,j,k), u_arr(i,j,k), u_arr(i+1,j,k), u_arr(i+2,j,k),
                dphi, w_L, w_R, dwdphi);
            
            double f_face = w_R;
            double d2f_phi, d4f_phi, d2dphi2; // We only care about d2dphi2 right now!
            weno_ao_43_boundary(point_f(i-1), point_f(i), point_f(i+1), point_f(i+2), dphi, d2f_phi, d4f_phi, d2dphi2);
                            
            flx(i,j,k,0) = alpha_v * f_face - (dphi*dphi/24.0) * d2dphi2;
        });
    }
}

void get_rhs_azi(const amrex::MultiFab& u, amrex::MultiFab& Rhs, const amrex::Geometry& geom) {
    amrex::MultiFab flux_act(u.boxArray(), u.DistributionMap(), 1, 1); 
    compute_flux_azi(u, flux_act, geom);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& f_act = flux_act.const_array(mfi);
        auto const& rhs = Rhs.array(mfi);
        
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double div_F = (f_act(i,j,k,0) - f_act(i-1,j,k,0)) / dphi;
            double factor = (m == 1) ? (1.0 / r_fixed) : (1.0 / (r_fixed * std::sin(theta_fixed)));
            rhs(i,j,k) = -factor * div_F; 
        });
    }
}

void run_azimuthal_advection_3() {
    amrex::Print() << "Initializing AMReX Advection Solver (Azimuthal)...\n";
    
    amrex::Box domain(amrex::IntVect(AMREX_D_DECL(0,0,0)), amrex::IntVect(AMREX_D_DECL(N-1,0,0)));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, -1.0, -1.0)}, {AMREX_D_DECL(2.0*M_PI, 1.0, 1.0)});
    
    amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 0, 0)};
    amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

    amrex::BoxArray ba(domain);
    ba.maxSize(16);
    amrex::DistributionMapping dm(ba);

    amrex::MultiFab u(ba, dm, 1, 3);
    amrex::MultiFab u1(ba, dm, 1, 3);
    amrex::MultiFab u2(ba, dm, 1, 3);
    amrex::MultiFab rhs(ba, dm, 1, 0);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& arr = u.array(mfi);
        amrex::ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            arr(i,j,k) = exact_sol_azi(get_phi(i), 0.0);
        });
    }
    apply_boundaries_azi(u, geom); 

    double t = 0.0, t_end = 1.0, CFL = 0.4;
    double factor = (m == 1) ? (1.0 / r_fixed) : (1.0 / (r_fixed * std::sin(theta_fixed)));
    double c_speed = alpha_v * factor;

    amrex::Print() << "Starting Time Integration...\n";
    BL_PROFILE_VAR("Evolution_Loop", pmain);
    int step = 0;
    double dt = CFL * dphi / c_speed;

    while (t < t_end) {
        if (t + dt > t_end) dt = t_end - t;

        get_rhs_azi(u, rhs, geom);
        amrex::MultiFab::Copy(u1, u, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u1, dt, rhs, 0, 0, 1, 0);
        apply_boundaries_azi(u1, geom);

        get_rhs_azi(u1, rhs, geom);
        amrex::MultiFab::LinComb(u2, 0.75, u, 0, 0.25, u1, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u2, 0.25 * dt, rhs, 0, 0, 1, 0);
        apply_boundaries_azi(u2, geom);

        get_rhs_azi(u2, rhs, geom);
        amrex::MultiFab::LinComb(u, 1.0/3.0, u, 0, 2.0/3.0, u2, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u, (2.0/3.0) * dt, rhs, 0, 0, 1, 0);
        apply_boundaries_azi(u, geom);
        
        t += dt; step++;
        if (step % 10 == 0) amrex::Print() << "Step " << step << " | Time: " << t << " / " << t_end << "\n";
    }
    BL_PROFILE_VAR_STOP(pmain); 

    
    amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum> reduce_ops;
    amrex::ReduceData<double, double, double, double, double> reduce_data(reduce_ops);
    using ReduceTuple = typename amrex::ReduceData<double, double, double, double, double>::Type;

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
    auto const& arr = u.const_array(mfi);
    reduce_ops.eval(mfi.validbox(), reduce_data,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
            double exact = exact_sol_azi(get_phi(i), t);
            double err   = std::abs(arr(i,j,k) - exact);

            double r_c = get_r(i);
            double w   = std::pow(r_c, (double)m) * dr;   // r^2 dr for m = 2

            return {err, err * w, err * err * w, err, w};
        });
    }

    ReduceTuple hv = reduce_data.value();
    double linf   = amrex::get<0>(hv);
    double num    = amrex::get<1>(hv);
    double l2num  = amrex::get<2>(hv);
    double sum    = amrex::get<3>(hv);
    double wsum   = amrex::get<4>(hv);

    amrex::ParallelDescriptor::ReduceRealMax(linf);
    amrex::ParallelDescriptor::ReduceRealSum(num);
    amrex::ParallelDescriptor::ReduceRealSum(l2num);
    amrex::ParallelDescriptor::ReduceRealSum(sum);
    amrex::ParallelDescriptor::ReduceRealSum(wsum);

    double l1_weighted = num / wsum;              // L1 variant for AFD
    double l1_plain    = sum / (double)N;    // unweighted
    double l2_weighted = std::sqrt(l2num / wsum);
    amrex::Print() << "\n--- Advection Benchmark Verification ---\n" << std::scientific << std::setprecision(8)
    << "L1 Error:    " << l1_plain << "\nL1 Error (Unweighted):    " << l1_plain << "\nL-inf Error: " << linf << "\nN: " << N << "\n";
    
    amrex::WriteSingleLevelPlotfile("plt_final", u, {"u"}, geom, t, step);
}

// ==============================================================================
// FIFTH ORDER
// ==============================================================================

void compute_flux_azi_5(const amrex::MultiFab& u, amrex::MultiFab& flux, const amrex::Geometry& geom) {
    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        amrex::Box flux_box = mfi.validbox(); flux_box.growLo(0, 1); 

        auto const& u_arr = u.const_array(mfi); 
        auto const& flx = flux.array(mfi);

        amrex::ParallelFor(flux_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double w_L, w_R, dwdphi;

            auto point_f = [&](int idx) { 
                return alpha_v * u_arr(idx,j,k); 
            };
            
            weno_ao_5_3_interpolation(
                u_arr(i-2,j,k), u_arr(i-1,j,k), u_arr(i,j,k), u_arr(i+1,j,k), u_arr(i+2,j,k),
                dphi, w_L, w_R, dwdphi);
            
            double f_face = w_R;
            
            double d2dphi2, d4dphi4; 
            weno_ao_63_boundary(
                point_f(i-2), point_f(i-1), point_f(i), 
                point_f(i+1), point_f(i+2), point_f(i+3), 
                dphi, d2dphi2, d4dphi4);

            flx(i,j,k,0) = alpha_v * f_face 
                         - (dphi * dphi / 24.0) * d2dphi2
                         + (7.0 * dphi * dphi * dphi * dphi / 5760.0) * d4dphi4;
        });
    }
}

void get_rhs_azi_5(const amrex::MultiFab& u, amrex::MultiFab& Rhs, const amrex::Geometry& geom) {
    amrex::MultiFab flux_act(u.boxArray(), u.DistributionMap(), 1, 1); 
    compute_flux_azi_5(u, flux_act, geom);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& f_act = flux_act.const_array(mfi);
        auto const& rhs = Rhs.array(mfi);
        
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            double div_F = (f_act(i,j,k,0) - f_act(i-1,j,k,0)) / dphi;
            double factor = (m == 1) ? (1.0 / r_fixed) : (1.0 / (r_fixed * std::sin(theta_fixed)));
            rhs(i,j,k) = -factor * div_F; 
        });
    }
}

void run_azimuthal_advection_5() {
    amrex::Print() << "Initializing AMReX Advection Solver (Azimuthal 5th Order)...\n";
    
    amrex::Box domain(amrex::IntVect(AMREX_D_DECL(0,0,0)), amrex::IntVect(AMREX_D_DECL(N-1,0,0)));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, -1.0, -1.0)}, {AMREX_D_DECL(2.0*M_PI, 1.0, 1.0)});
    
    amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 0, 0)};
    amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

    amrex::BoxArray ba(domain);
    ba.maxSize(16);
    amrex::DistributionMapping dm(ba);

    amrex::MultiFab u (ba, dm, 1, 3);
    amrex::MultiFab u0(ba, dm, 1, 3);
    amrex::MultiFab u1(ba, dm, 1, 3);
    amrex::MultiFab u2(ba, dm, 1, 3);
    amrex::MultiFab u3(ba, dm, 1, 3);
    amrex::MultiFab u4(ba, dm, 1, 3);
    amrex::MultiFab rhs(ba, dm, 1, 0);

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
        auto const& arr = u.array(mfi);
        amrex::ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            arr(i,j,k) = exact_sol_azi(get_phi(i), 0.0);
        });
    }
    apply_boundaries_azi(u, geom); 

    double t = 0.0, t_end = 1.0;
    double factor = (m == 1) ? (1.0 / r_fixed) : (1.0 / (r_fixed * std::sin(theta_fixed)));
    double c_speed = alpha_v * factor;

    double N_base = 64.0;
    double base_CFL = 0.4;
    double scaled_CFL = base_CFL * std::pow(N_base / (double)N, 0.25);
    double dt = scaled_CFL * dphi / c_speed;

    amrex::Print() << "Starting Time Integration... Scaled CFL = " << scaled_CFL << "\n";
    BL_PROFILE_VAR("Evolution_Loop", pmain);
    int step = 0;

    while (t < t_end) {
        if (t + dt > t_end) dt = t_end - t;

        amrex::MultiFab::Copy(u0, u, 0, 0, 1, 3);

        // --- Stage 1 ---
        get_rhs_azi_5(u0, rhs, geom);
        amrex::MultiFab::Copy(u1, u0, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u1, dt * 0.39175222700392, rhs, 0, 0, 1, 0);
        apply_boundaries_azi(u1, geom);

        // --- Stage 2 ---
        get_rhs_azi_5(u1, rhs, geom);
        amrex::MultiFab::LinComb(u2, 0.44437049406734, u0, 0, 0.55562950593266, u1, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u2, dt * 0.36841059262959, rhs, 0, 0, 1, 0);
        apply_boundaries_azi(u2, geom);

        // --- Stage 3 ---
        get_rhs_azi_5(u2, rhs, geom);
        amrex::MultiFab::LinComb(u3, 0.62010185138540, u0, 0, 0.37989814861460, u2, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u3, dt * 0.25189177424738, rhs, 0, 0, 1, 0);
        apply_boundaries_azi(u3, geom);

        // --- Stage 4 ---
        amrex::MultiFab rhs_u3(ba, dm, 1, 0); 
        get_rhs_azi_5(u3, rhs_u3, geom);
        amrex::MultiFab::LinComb(u4, 0.17807995410773, u0, 0, 0.82192004589227, u3, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u4, dt * 0.54497475021237, rhs_u3, 0, 0, 1, 0);
        apply_boundaries_azi(u4, geom);

        // --- Stage 5 ---
        get_rhs_azi_5(u4, rhs, geom);
        u.setVal(0.0);
        amrex::MultiFab::Saxpy(u, 0.00683325884039, u0, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u, 0.51723167208978, u2, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u, 0.12759831133288, u3, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u, 0.34833675773694, u4, 0, 0, 1, 0);

        amrex::MultiFab::Saxpy(u, dt * 0.08460416338212, rhs_u3, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(u, dt * 0.22600748319395, rhs, 0, 0, 1, 0);   
        apply_boundaries_azi(u, geom);

        t += dt; 
        step++;
        if (step % 10 == 0) amrex::Print() << "Step " << step << " | Time: " << t << " / " << t_end << "\n";
    }
    BL_PROFILE_VAR_STOP(pmain); 

    
    amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum> reduce_ops;
    amrex::ReduceData<double, double, double, double, double> reduce_data(reduce_ops);
    using ReduceTuple = typename amrex::ReduceData<double, double, double, double, double>::Type;

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
    auto const& arr = u.const_array(mfi);
    reduce_ops.eval(mfi.validbox(), reduce_data,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
            double exact = exact_sol_azi(get_phi(i), t);
            double err   = std::abs(arr(i,j,k) - exact);

            double r_c = get_r(i);
            double w   = std::pow(r_c, (double)m) * dr;   // r^2 dr for m = 2

            return {err, err * w, err * err * w, err, w};
        });
    }

    ReduceTuple hv = reduce_data.value();
    double linf   = amrex::get<0>(hv);
    double num    = amrex::get<1>(hv);
    double l2num  = amrex::get<2>(hv);
    double sum    = amrex::get<3>(hv);
    double wsum   = amrex::get<4>(hv);

    amrex::ParallelDescriptor::ReduceRealMax(linf);
    amrex::ParallelDescriptor::ReduceRealSum(num);
    amrex::ParallelDescriptor::ReduceRealSum(l2num);
    amrex::ParallelDescriptor::ReduceRealSum(sum);
    amrex::ParallelDescriptor::ReduceRealSum(wsum);

    double l1_weighted = num / wsum;              // L1 variant for AFD
    double l1_plain    = sum / (double)N;    // unweighted
    double l2_weighted = std::sqrt(l2num / wsum);
    amrex::Print() << "\n--- Advection Benchmark Verification ---\n" << std::scientific << std::setprecision(8)
    << "L1 Error:    " << l1_plain << "\nL1 Error (Unweighted):    " << l1_plain << "\nL-inf Error: " << linf << "\nN: " << N << "\n";
    
    amrex::WriteSingleLevelPlotfile("plt_final", u, {"u"}, geom, t, step);
}

// ==============================================================================
// MAIN RUNNER
// ==============================================================================

int main(int argc, char* argv[]) {
    amrex::Initialize(argc, argv);
    
    if constexpr (ORDER == 3) {
        if constexpr (TEST_MODE == Mode::RADIAL)           run_radial_advection_3();
        else if constexpr (TEST_MODE == Mode::MERIDIONAL)  run_meridional_advection_3();
        else if constexpr (TEST_MODE == Mode::Z)           run_z_advection_3();
        else if constexpr (TEST_MODE == Mode::AZIMUTHAL)   run_azimuthal_advection_3();
    } 

    else if constexpr (ORDER == 5) {
        if constexpr (TEST_MODE == Mode::RADIAL)           run_radial_advection_5();
        else if constexpr (TEST_MODE == Mode::MERIDIONAL)  run_meridional_advection_5();
        else if constexpr (TEST_MODE == Mode::Z)           run_z_advection_5();
        else if constexpr (TEST_MODE == Mode::AZIMUTHAL)   run_azimuthal_advection_5();
    } 

    amrex::Finalize();
    return 0;
}
