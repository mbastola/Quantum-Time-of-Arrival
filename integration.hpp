#ifndef QTOA_INTEGRATION_HPP
#define QTOA_INTEGRATION_HPP

#include "common.hpp"

namespace qtoa {

template <typename F>
double integrate_real(F f, double a, double b, int steps) {
    if (steps < 2) {
        steps = 2;
    }
    if (steps % 2 != 0) {
        ++steps;
    }

    const double h = (b - a) / static_cast<double>(steps);
    double sum = f(a) + f(b);
    for (int i = 1; i < steps; ++i) {
        const double x = a + static_cast<double>(i) * h;
        sum += (i % 2 == 0 ? 2.0 : 4.0) * f(x);
    }
    return sum * h / 3.0;
}

template <typename F>
Complex integrate_complex(F f, double a, double b, int steps) {
    if (steps < 2) {
        steps = 2;
    }
    if (steps % 2 != 0) {
        ++steps;
    }

    const double h = (b - a) / static_cast<double>(steps);
    Complex sum = f(a) + f(b);
    for (int i = 1; i < steps; ++i) {
        const double x = a + static_cast<double>(i) * h;
        sum += (i % 2 == 0 ? 2.0 : 4.0) * f(x);
    }
    return sum * h / 3.0;
}

}  // namespace qtoa

#endif
