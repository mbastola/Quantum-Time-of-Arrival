#include "physics.hpp"

#include <cmath>
#include <stdexcept>

namespace qtoa {

Complex psi0(const Params& p, double x, double t) {
    const double d2 = p.delta * p.delta;
    const double prefactor = std::pow(d2 / (2.0 * PI), 0.25) *
                             std::exp(-p.k0 * p.k0 * d2);
    const Complex spreading = d2 + I * t / 2.0;
    const Complex numerator = 2.0 * d2 * p.k0 + I * (x - p.x0);
    const Complex denominator = 4.0 * d2 + 2.0 * I * t;
    return prefactor / std::sqrt(spreading) *
           std::exp((numerator * numerator) / denominator);
}

double psi0_density(const Params& p, double x, double t) {
    return std::norm(psi0(p, x, t));
}

Complex propagator(double xf, double tf, double xi, double ti) {
    const double dt = tf - ti;
    if (std::fabs(dt) < 1.0e-14) {
        throw std::runtime_error("propagator is singular for equal times");
    }

    return 1.0 / std::sqrt(2.0 * PI * I * dt) *
           std::exp(I * (xf - xi) * (xf - xi) / (2.0 * dt));
}

}  // namespace qtoa
