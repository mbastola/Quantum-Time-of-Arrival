#ifndef QTOA_PHYSICS_HPP
#define QTOA_PHYSICS_HPP

#include "common.hpp"

namespace qtoa {

Complex psi0(const Params& p, double x, double t);
double psi0_density(const Params& p, double x, double t);
Complex propagator(double xf, double tf, double xi, double ti);

}  // namespace qtoa

#endif
