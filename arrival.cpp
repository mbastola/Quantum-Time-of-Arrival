#include "arrival.hpp"

#include "integration.hpp"
#include "physics.hpp"

#include <algorithm>
#include <cmath>
#include <random>

namespace qtoa {
namespace {

double random_x(const Params& p, double detector_width, std::mt19937& rng) {
    std::uniform_real_distribution<double> dist(
        p.xd - detector_width / 2.0, p.xd + detector_width / 2.0);
    return dist(rng);
}

void evaluate_dynamic_sequence(const Params& p,
                               const std::vector<double>& seq,
                               double detector_width,
                               double& numerator,
                               double& denominator) {
    const int n = static_cast<int>(seq.size());
    std::vector<double> a(n + 1, 0.0);
    std::vector<double> b(n + 1, 1.0);
    std::vector<std::vector<Complex> > alpha(
        n + 1, std::vector<Complex>(n + 1, Complex(0.0, 0.0)));

    for (int j = 1; j <= n; ++j) {
        alpha[0][j] = psi0(p, p.xd, seq[j - 1]);
    }

    a[1] = std::abs(alpha[0][1]) * std::sqrt(detector_width);
    b[1] = safe_sqrt(1.0 - a[1] * a[1]);

    for (int i = 1; i <= n - 1; ++i) {
        for (int j = i + 1; j <= n; ++j) {
            if (std::fabs(seq[i - 1] - seq[j - 1]) < 1.0e-14) {
                alpha[i][j] = Complex(0.0, 0.0);
            } else if (b[i] == 0.0) {
                alpha[i][j] = Complex(0.0, 0.0);
            } else {
                const Complex beta =
                    propagator(p.xd, 0.0, p.xd, seq[j - 1] - seq[i - 1]);
                alpha[i][j] =
                    (alpha[i - 1][j] -
                     alpha[i - 1][i] * beta * detector_width) /
                    b[i];
            }
        }
        a[i + 1] = std::abs(alpha[i][i + 1]) *
                   std::sqrt(detector_width);
        b[i + 1] = safe_sqrt(1.0 - a[i + 1] * a[i + 1]);
    }

    double survival = 1.0;
    for (int j = 1; j <= n; ++j) {
        const double prob = survival * a[j] * a[j];
        numerator += seq[j - 1] * prob;
        denominator += prob;
        survival *= b[j] * b[j];
    }
}

void recurse_sequences(const Params& p,
                       const std::vector<double>& grid,
                       int measurements,
                       int depth,
                       bool directional,
                       int start,
                       double detector_width,
                       std::vector<double>& seq,
                       double& numerator,
                       double& denominator) {
    if (depth == measurements) {
        evaluate_dynamic_sequence(p, seq, detector_width, numerator,
                                  denominator);
        return;
    }

    const int first = directional ? start : 0;
    for (int i = first; i < static_cast<int>(grid.size()); ++i) {
        seq[depth] = grid[i];
        recurse_sequences(p, grid, measurements, depth + 1, directional, i,
                          detector_width, seq, numerator, denominator);
    }
}

}  // namespace

double exact_a1(const Params& p, double detector_width, int steps) {
    const double lo = p.xd - detector_width / 2.0;
    const double hi = p.xd + detector_width / 2.0;
    const double prob = integrate_real(
        [&](double x) { return psi0_density(p, x, p.ti); }, lo, hi, steps);
    return safe_sqrt(prob);
}

double small_a1(const Params& p, double detector_width) {
    return std::abs(psi0(p, p.xd, p.ti)) * std::sqrt(detector_width);
}

double exact_a2(const Params& p, double detector_width, int n_for_times) {
    Params local = p;
    local.n = n_for_times;
    const std::vector<double> t = make_times(local.n, local.ti, local.tf);
    const double t1 = t[1];
    const double t2 = t[2];
    const double lo = local.xd - detector_width / 2.0;
    const double hi = local.xd + detector_width / 2.0;

    const double c1 = exact_a1(local, detector_width, 160);
    const double b1 = safe_sqrt(1.0 - c1 * c1);
    if (c1 == 0.0 || b1 == 0.0) {
        return 0.0;
    }

    const int inner_steps = 120;
    const int outer_steps = 120;
    const double prob = integrate_real(
        [&](double x) {
            const Complex phi1 = integrate_complex(
                                     [&](double xp) {
                                         return psi0(local, xp, t1) *
                                                propagator(x, t2, xp, t1);
                                     },
                                     lo, hi, inner_steps) /
                                 c1;
            const Complex psi1 = (psi0(local, x, t2) - c1 * phi1) / b1;
            return std::norm(psi1);
        },
        lo, hi, outer_steps);

    return safe_sqrt(prob);
}

ArrivalResult small_detector(const Params& p, int n, double detector_width) {
    ArrivalResult r;
    r.n = n;
    r.delta1 = detector_width;
    r.tbar = 0.0;
    r.t = make_times(n, p.ti, p.tf);
    r.a.assign(n + 1, 0.0);
    r.b.assign(n + 1, 1.0);
    r.prob.assign(n + 1, 0.0);
    r.alpha.assign(n + 1, std::vector<Complex>(n + 2, Complex(0.0, 0.0)));

    std::vector<std::vector<Complex> > beta(
        n + 1, std::vector<Complex>(n + 2, Complex(0.0, 0.0)));

    for (int j = 1; j <= n + 1; ++j) {
        r.alpha[0][j] = psi0(p, p.xd, r.t[j]);
    }

    for (int i = 1; i <= n; ++i) {
        for (int j = i + 1; j <= n + 1; ++j) {
            beta[i][j] = propagator(p.xd, 0.0, p.xd, r.t[j] - r.t[i]);
        }
    }

    r.a[1] = std::abs(r.alpha[0][1]) * std::sqrt(detector_width);
    r.b[1] = safe_sqrt(1.0 - r.a[1] * r.a[1]);

    for (int i = 1; i <= n - 1; ++i) {
        for (int j = i + 1; j <= n + 1; ++j) {
            if (r.b[i] == 0.0) {
                r.alpha[i][j] = Complex(0.0, 0.0);
            } else {
                r.alpha[i][j] =
                    (r.alpha[i - 1][j] -
                     r.alpha[i - 1][i] * beta[i][j] * detector_width) /
                    r.b[i];
            }
        }
        r.a[i + 1] = std::abs(r.alpha[i][i + 1]) *
                     std::sqrt(detector_width);
        r.b[i + 1] = safe_sqrt(1.0 - r.a[i + 1] * r.a[i + 1]);
    }

    double survival = 1.0;
    for (int i = 1; i <= n; ++i) {
        r.prob[i] = survival * r.a[i] * r.a[i];
        survival *= r.b[i] * r.b[i];
    }

    double top = 0.0;
    double bottom = 0.0;
    for (int i = 1; i <= n; ++i) {
        top += r.t[i] * r.prob[i];
        bottom += r.prob[i];
    }
    r.tbar = bottom == 0.0 ? 0.0 : top / bottom;
    return r;
}

ArrivalResult random_detector(const Params& p,
                              int n,
                              double detector_width,
                              int samples,
                              unsigned seed) {
    ArrivalResult r;
    r.n = n;
    r.delta1 = detector_width;
    r.tbar = 0.0;
    r.t = make_times(n, p.ti, p.tf);
    r.a.assign(n + 1, 0.0);
    r.b.assign(n + 1, 0.0);
    r.prob.assign(n + 1, 0.0);

    std::mt19937 rng(seed);

    for (int sample = 0; sample < samples; ++sample) {
        std::vector<double> a(n + 1, 0.0);
        std::vector<double> b(n + 1, 1.0);
        std::vector<std::vector<Complex> > alpha(
            n + 1, std::vector<Complex>(n + 2, Complex(0.0, 0.0)));
        std::vector<std::vector<Complex> > beta(
            n + 1, std::vector<Complex>(n + 2, Complex(0.0, 0.0)));

        a[1] = std::abs(psi0(p, random_x(p, detector_width, rng), r.t[1])) *
               std::sqrt(detector_width);
        b[1] = safe_sqrt(1.0 - a[1] * a[1]);

        for (int j = 1; j <= n + 1; ++j) {
            alpha[0][j] =
                psi0(p, random_x(p, detector_width, rng), r.t[j]);
        }

        for (int i = 1; i <= n; ++i) {
            for (int j = i + 1; j <= n + 1; ++j) {
                beta[i][j] =
                    propagator(random_x(p, detector_width, rng), 0.0,
                               random_x(p, detector_width, rng),
                               r.t[j] - r.t[i]);
            }
        }

        for (int i = 1; i <= n - 1; ++i) {
            for (int j = i + 1; j <= n + 1; ++j) {
                if (b[i] == 0.0) {
                    alpha[i][j] = Complex(0.0, 0.0);
                } else {
                    alpha[i][j] =
                        (alpha[i - 1][j] -
                         alpha[i - 1][i] * beta[i][j] * detector_width) /
                        b[i];
                }
            }
            a[i + 1] = std::abs(alpha[i][i + 1]) *
                       std::sqrt(detector_width);
            b[i + 1] = safe_sqrt(1.0 - a[i + 1] * a[i + 1]);
        }

        for (int i = 1; i <= n; ++i) {
            r.a[i] += a[i];
            r.b[i] += b[i];
        }
    }

    for (int i = 1; i <= n; ++i) {
        r.a[i] /= static_cast<double>(samples);
        r.b[i] /= static_cast<double>(samples);
    }

    double survival = 1.0;
    for (int i = 1; i <= n; ++i) {
        r.prob[i] = survival * r.a[i] * r.a[i];
        survival *= r.b[i] * r.b[i];
    }

    double top = 0.0;
    double bottom = 0.0;
    for (int i = 1; i <= n; ++i) {
        top += r.t[i] * r.prob[i];
        bottom += r.prob[i];
    }
    r.tbar = bottom == 0.0 ? 0.0 : top / bottom;
    return r;
}

Complex evolved_state(const Params& p,
                      const ArrivalResult& r,
                      double x,
                      double time,
                      int state_index) {
    state_index = std::max(0, std::min(state_index, r.n));
    Complex value = psi0(p, x, time);

    for (int i = 1; i <= state_index; ++i) {
        if (r.b[i] == 0.0) {
            return Complex(0.0, 0.0);
        }
        const Complex before_measurement = r.alpha[i - 1][i];
        value = (value -
                 before_measurement *
                     propagator(x, 0.0, p.xd, time - r.t[i]) * r.delta1) /
                r.b[i];
    }
    return value;
}

int state_index_for_time(const ArrivalResult& r, double time) {
    if (r.n < 2) {
        return 0;
    }
    const double dt = r.t[2] - r.t[1];
    int index = static_cast<int>(((time - r.t[1]) + dt) / dt);
    return std::max(0, std::min(index, r.n));
}

double dynamic_arrival_time(const Params& p,
                            int time_steps,
                            int measurements,
                            double detector_width,
                            bool directional) {
    const std::vector<double> grid = make_time_grid(time_steps, p.ti, p.tf);
    std::vector<double> seq(measurements, 0.0);
    double numerator = 0.0;
    double denominator = 0.0;

    recurse_sequences(p, grid, measurements, 0, directional, 0, detector_width,
                      seq, numerator, denominator);

    return denominator == 0.0 ? 0.0 : numerator / denominator;
}

}  // namespace qtoa
