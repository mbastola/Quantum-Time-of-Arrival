#ifndef QTOA_ARRIVAL_HPP
#define QTOA_ARRIVAL_HPP

#include "common.hpp"

namespace qtoa {

double exact_a1(const Params& p, double detector_width, int steps);
double small_a1(const Params& p, double detector_width);
double exact_a2(const Params& p, double detector_width, int n_for_times);

ArrivalResult small_detector(const Params& p, int n, double detector_width);
ArrivalResult random_detector(const Params& p,
                              int n,
                              double detector_width,
                              int samples,
                              unsigned seed);

Complex evolved_state(const Params& p,
                      const ArrivalResult& r,
                      double x,
                      double time,
                      int state_index);
int state_index_for_time(const ArrivalResult& r, double time);

double dynamic_arrival_time(const Params& p,
                            int time_steps,
                            int measurements,
                            double detector_width,
                            bool directional);

}  // namespace qtoa

#endif
