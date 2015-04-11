#ifndef QTOA_COMMON_HPP
#define QTOA_COMMON_HPP

#include <complex>
#include <string>
#include <vector>

namespace qtoa {

typedef std::complex<double> Complex;

extern const double PI;
extern const Complex I;

struct Params {
    double k0;
    double x0;
    double xd;
    double delta;
    double delta1;
    double ti;
    double tf;
    int n;

    Params()
        : k0(20.0),
          x0(-5.0),
          xd(-3.0),
          delta(0.5),
          delta1(0.005),
          ti(-0.02),
          tf(0.22),
          n(100) {}
};

struct Point {
    double x;
    double y;
};

struct Series {
    std::string name;
    std::string color;
    bool points;
    std::vector<Point> values;
};

struct ArrivalResult {
    int n;
    double delta1;
    double tbar;
    std::vector<double> t;
    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> prob;
    std::vector<std::vector<Complex> > alpha;
};

double clamp_nonnegative(double x);
double safe_sqrt(double x);
std::vector<double> make_times(int n, double ti, double tf);
std::vector<double> make_time_grid(int count, double ti, double tf);

}  // namespace qtoa

#endif
