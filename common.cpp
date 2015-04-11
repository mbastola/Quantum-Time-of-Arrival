#include "common.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace qtoa {

extern const double PI = std::acos(-1.0);
extern const Complex I = Complex(0.0, 1.0);

double clamp_nonnegative(double x) {
    if (x < 0.0 && x > -1.0e-12) {
        return 0.0;
    }
    return x;
}

double safe_sqrt(double x) {
    return std::sqrt(std::max(0.0, clamp_nonnegative(x)));
}

std::vector<double> make_times(int n, double ti, double tf) {
    if (n < 2) {
        throw std::runtime_error("n must be at least 2");
    }

    std::vector<double> t(n + 2, 0.0);
    const double dt = (tf - ti) / static_cast<double>(n - 1);
    t[0] = ti - dt;
    for (int i = 1; i <= n; ++i) {
        t[i] = ti + static_cast<double>(i - 1) * dt;
    }
    t[n + 1] = tf + dt;
    return t;
}

std::vector<double> make_time_grid(int count, double ti, double tf) {
    if (count < 2) {
        throw std::runtime_error("time grid count must be at least 2");
    }

    std::vector<double> t(count);
    const double dt = (tf - ti) / static_cast<double>(count - 1);
    for (int i = 0; i < count; ++i) {
        t[i] = ti + static_cast<double>(i) * dt;
    }
    return t;
}

}  // namespace qtoa
