#include "generation.hpp"

#include "arrival.hpp"
#include "physics.hpp"
#include "plotting.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

namespace qtoa {
namespace {

double random_a1(const Params& p,
                 double detector_width,
                 int samples,
                 unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(
        p.xd - detector_width / 2.0, p.xd + detector_width / 2.0);

    double total = 0.0;
    for (int i = 0; i < samples; ++i) {
        total += std::abs(psi0(p, dist(rng), p.ti)) *
                 std::sqrt(detector_width);
    }
    return total / static_cast<double>(samples);
}

void generate_a1_plots(const Params& p, const std::string& out_dir) {
    std::vector<Point> exact;
    std::vector<Point> small;
    std::vector<Point> random;
    std::vector<Point> diff_small;
    std::vector<Point> diff_random;

    int index = 0;
    for (double d = 0.00005; d <= 0.5000001; d += 0.005) {
        const double a_exact = exact_a1(p, d, 240);
        const double a_small = small_a1(p, d);
        const double a_random = random_a1(p, d, 1200, 2000u + index);
        exact.push_back(Point{d, a_exact});
        small.push_back(Point{d, a_small});
        random.push_back(Point{d, a_random});
        diff_small.push_back(Point{d, std::fabs(a_exact - a_small)});
        diff_random.push_back(Point{d, std::fabs(a_exact - a_random)});
        ++index;
    }

    write_plot(out_dir, "fig3_1_a1_comparison",
               "Fig. 3.1: a1 by detector width", "delta1", "a1",
               std::vector<Series>{
                   make_series("exact integral", "#1f4e79", false, exact),
                   make_series("small detector", "#c43c39", false, small),
                   make_series("random detector", "#2a9d55", true, random)});

    write_plot(out_dir, "fig3_2_a1_difference",
               "Fig. 3.2: difference from exact a1", "delta1",
               "absolute difference",
               std::vector<Series>{
                   make_series("|exact-small|", "#c43c39", false, diff_small),
                   make_series("|exact-random|", "#2a9d55", true,
                               diff_random)});
}

void generate_a2_plot(const Params& p, const std::string& out_dir) {
    std::vector<Point> exact;
    std::vector<Point> small;
    std::vector<Point> random;

    int index = 0;
    for (double d = 0.00005; d <= 0.5000001; d += 0.01) {
        const double a_exact = exact_a2(p, d, 100);
        const double a_small = small_detector(p, 100, d).a[2];
        const double a_random =
            random_detector(p, 2, d, 1000, 3000u + index).a[2];
        exact.push_back(Point{d, a_exact});
        small.push_back(Point{d, a_small});
        random.push_back(Point{d, a_random});
        ++index;
    }

    write_plot(out_dir, "fig3_3_a2_comparison",
               "Fig. 3.3: a2 by detector width", "delta1", "a2",
               std::vector<Series>{
                   make_series("numeric integral", "#1f4e79", false, exact),
                   make_series("small detector", "#c43c39", false, small),
                   make_series("random detector", "#2a9d55", true, random)});
}

void generate_an_probability_plots(const Params& p,
                                   const std::string& out_dir,
                                   bool quick) {
    const int n = 100;
    {
        const double d = 0.005;
        const ArrivalResult small_result = small_detector(p, n, d);
        const ArrivalResult random_result =
            random_detector(p, n, d, quick ? 40 : 100, 4000u);

        write_plot(out_dir, "fig3_4_an_small_width",
                   "Fig. 3.4: amplitudes for delta1=0.005", "t_n", "a_n",
                   std::vector<Series>{
                       make_series("small detector", "#c43c39", false,
                                   points_from_arrival(small_result,
                                                       small_result.a)),
                       make_series("random detector", "#2a9d55", false,
                                   points_from_arrival(random_result,
                                                       random_result.a))});

        write_plot(out_dir, "fig3_5_probability_small_width",
                   "Fig. 3.5: detection probability for delta1=0.005", "t_n",
                   "pi_n",
                   std::vector<Series>{
                       make_series("small detector", "#c43c39", false,
                                   points_from_arrival(small_result,
                                                       small_result.prob)),
                       make_series("random detector", "#2a9d55", false,
                                   points_from_arrival(random_result,
                                                       random_result.prob))});
    }

    {
        const double d = 0.035;
        const ArrivalResult small_result = small_detector(p, n, d);
        const ArrivalResult random_result =
            random_detector(p, n, d, quick ? 200 : 2000, 5000u);

        write_plot(out_dir, "fig3_6_an_large_width",
                   "Fig. 3.6: amplitudes for delta1=0.035", "t_n", "a_n",
                   std::vector<Series>{
                       make_series("small detector", "#c43c39", false,
                                   points_from_arrival(small_result,
                                                       small_result.a)),
                       make_series("random detector", "#2a9d55", false,
                                   points_from_arrival(random_result,
                                                       random_result.a))});

        write_plot(out_dir, "fig3_7_probability_large_width",
                   "Fig. 3.7: detection probability for delta1=0.035", "t_n",
                   "pi_n",
                   std::vector<Series>{
                       make_series("small detector", "#c43c39", false,
                                   points_from_arrival(small_result,
                                                       small_result.prob)),
                       make_series("random detector", "#2a9d55", false,
                                   points_from_arrival(random_result,
                                                       random_result.prob))});
    }
}

void generate_tbar_plots(const Params& p,
                         const std::string& out_dir,
                         bool quick) {
    std::vector<Point> small_width;
    std::vector<Point> random_width;
    std::vector<Point> classical;

    int index = 0;
    for (double d = 0.00005; d <= 0.0500001; d += 0.001) {
        const ArrivalResult small_result = small_detector(p, 100, d);
        const ArrivalResult random_result =
            random_detector(p, 100, d, quick ? 4 : 10, 6000u + index);
        small_width.push_back(Point{d, small_result.tbar});
        random_width.push_back(Point{d, random_result.tbar});
        classical.push_back(Point{d, (p.xd - p.x0) / p.k0});
        ++index;
    }

    write_plot(out_dir, "fig3_8_tbar_by_delta1",
               "Fig. 3.8: arrival time by detector width", "delta1", "tbar",
               std::vector<Series>{
                   make_series("small detector", "#c43c39", false,
                               small_width),
                   make_series("random detector", "#2a9d55", false,
                               random_width),
                   make_series("classical t=0.1", "#333333", false,
                               classical)});

    std::vector<Point> n_d005;
    std::vector<Point> n_d05;
    std::vector<Point> n_d5;
    std::vector<Point> n_classical;
    for (int n = 3; n <= 100; ++n) {
        n_d005.push_back(Point{static_cast<double>(n),
                               small_detector(p, n, 0.005).tbar});
        n_d05.push_back(
            Point{static_cast<double>(n), small_detector(p, n, 0.05).tbar});
        n_d5.push_back(
            Point{static_cast<double>(n), small_detector(p, n, 0.5).tbar});
        n_classical.push_back(Point{static_cast<double>(n),
                                    (p.xd - p.x0) / p.k0});
    }

    write_plot(out_dir, "fig3_9_tbar_small_by_n",
               "Fig. 3.9: small-detector tbar by n", "n", "tbar",
               std::vector<Series>{
                   make_series("delta1=0.005", "#1f4e79", false, n_d005),
                   make_series("delta1=0.05", "#c43c39", false, n_d05),
                   make_series("delta1=0.5", "#2a9d55", false, n_d5),
                   make_series("classical t=0.1", "#333333", false,
                               n_classical)});

    std::vector<Point> random_d05;
    std::vector<Point> random_d5;
    std::vector<Point> random_classical;
    for (int n = 5; n <= 100; n += 5) {
        random_d05.push_back(Point{
            static_cast<double>(n),
            random_detector(p, n, 0.05, quick ? 20 : 100, 7000u + n).tbar});
        random_d5.push_back(Point{
            static_cast<double>(n),
            random_detector(p, n, 0.5, quick ? 80 : 400, 8000u + n).tbar});
        random_classical.push_back(Point{static_cast<double>(n),
                                         (p.xd - p.x0) / p.k0});
    }

    write_plot(out_dir, "fig3_10_tbar_random_by_n",
               "Fig. 3.10: random-detector tbar by n", "n", "tbar",
               std::vector<Series>{
                   make_series("delta1=0.05", "#c43c39", false, random_d05),
                   make_series("delta1=0.5", "#2a9d55", false, random_d5),
                   make_series("classical t=0.1", "#333333", false,
                               random_classical)});
}

void generate_wavefunction_plot(const Params& p, const std::string& out_dir) {
    Params local = p;
    local.n = 10;
    local.delta1 = 0.005;
    const ArrivalResult r = small_detector(local, local.n, local.delta1);
    const double times_to_plot[] = {0.02, 0.10, 0.21};
    const char* names[] = {"t=0.02", "t=0.10", "t=0.21"};
    const char* colors[] = {"#1f4e79", "#c43c39", "#2a9d55"};

    std::vector<Series> series;
    for (int k = 0; k < 3; ++k) {
        const double t = times_to_plot[k];
        const int state = state_index_for_time(r, t);
        std::vector<Point> pts;
        for (int i = 0; i <= 300; ++i) {
            const double x = -7.0 + 5.0 * static_cast<double>(i) / 300.0;
            const Complex psi = evolved_state(local, r, x, t, state);
            pts.push_back(Point{x, std::norm(psi)});
        }
        series.push_back(make_series(names[k], colors[k], false, pts));
    }

    write_plot(out_dir, "wavefunction_evolution",
               "Post-measurement wavefunction evolution", "x", "|psi_n|^2",
               series);
}

void write_summary(const Params& p,
                   const std::string& out_dir,
                   bool quick) {
    std::ofstream out((out_dir + "/verification.txt").c_str());
    out << std::setprecision(12);
    out << "Quantum time-of-arrival C++11 verification\n\n";
    out << "Default thesis parameters:\n";
    out << "k0=" << p.k0 << ", x0=" << p.x0 << ", xd=" << p.xd
        << ", delta=" << p.delta << ", ti=" << p.ti << ", tf=" << p.tf
        << "\n";
    out << "classical arrival time = " << (p.xd - p.x0) / p.k0 << "\n\n";

    const ArrivalResult app3 = small_detector(p, 100, 0.005);
    out << "Appendix 3 small-detector check, n=100, delta1=0.005:\n";
    out << "tbar = " << app3.tbar << " (appendix text reports 0.100794)\n\n";

    out << "Chapter 4 comparison check, n=7, delta1=0.005:\n";
    out << "small-detector tbar = " << small_detector(p, 7, 0.005).tbar
        << " (thesis text reports 0.100817)\n\n";

    out << "Additional non-dynamic checks generated from the same recurrence, n=20:\n";
    out << "delta1=0.005 small tbar = " << small_detector(p, 20, 0.005).tbar
        << "\n";
    out << "delta1=0.035 small tbar = " << small_detector(p, 20, 0.035).tbar
        << "\n";
    if (!quick) {
        out << "delta1=0.005 random tbar, m=100 = "
            << random_detector(p, 20, 0.005, 100, 9100u).tbar << "\n";
        out << "delta1=0.035 random tbar, m=10000 = "
            << random_detector(p, 20, 0.035, 10000, 9200u).tbar << "\n";
    } else {
        out << "random Table 3.1 checks skipped/reduced in --quick mode\n";
    }

    out << "\nTable 4.1 dynamic-time checks:\n";
    const struct {
        int m;
        int n;
    } rows[] = {{10, 4}, {7, 7}, {20, 4}, {100, 3}};
    for (std::size_t i = 0; i < sizeof(rows) / sizeof(rows[0]); ++i) {
        const double nd =
            dynamic_arrival_time(p, rows[i].m, rows[i].n, 0.005, false);
        out << "m=" << rows[i].m << ", n=" << rows[i].n
            << ", nondirectional tbar=" << nd;
        if (rows[i].n <= 4) {
            const double dir =
                dynamic_arrival_time(p, rows[i].m, rows[i].n, 0.005, true);
            out << ", directional tbar=" << dir;
            if (rows[i].m == 10 && rows[i].n == 4) {
                out << " (thesis table prints 0.100575; this recurrence gives 0.100475)";
            }
        }
        out << "\n";
    }
}

}  // namespace

void generate_all(const std::string& out_dir, bool quick) {
    Params p;
    ensure_dir(out_dir);

    std::cout << "Generating a1 plots...\n";
    generate_a1_plots(p, out_dir);
    std::cout << "Generating a2 plot...\n";
    generate_a2_plot(p, out_dir);
    std::cout << "Generating amplitude/probability plots...\n";
    generate_an_probability_plots(p, out_dir, quick);
    std::cout << "Generating tbar plots...\n";
    generate_tbar_plots(p, out_dir, quick);
    std::cout << "Generating wavefunction evolution plot...\n";
    generate_wavefunction_plot(p, out_dir);
    std::cout << "Writing verification summary...\n";
    write_summary(p, out_dir, quick);
}

}  // namespace qtoa
