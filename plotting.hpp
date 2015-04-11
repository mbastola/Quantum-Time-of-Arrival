#ifndef QTOA_PLOTTING_HPP
#define QTOA_PLOTTING_HPP

#include "common.hpp"

namespace qtoa {

void ensure_dir(const std::string& path);

Series make_series(const std::string& name,
                   const std::string& color,
                   bool points,
                   const std::vector<Point>& values);
std::vector<Point> points_from_arrival(const ArrivalResult& r,
                                       const std::vector<double>& y);

void write_plot(const std::string& out_dir,
                const std::string& stem,
                const std::string& title,
                const std::string& xlabel,
                const std::string& ylabel,
                const std::vector<Series>& series);

}  // namespace qtoa

#endif
