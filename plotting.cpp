#include "plotting.hpp"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <sys/stat.h>
#include <sys/types.h>

namespace qtoa {
namespace {

std::string escape_xml(const std::string& s) {
    std::string out;
    for (std::size_t i = 0; i < s.size(); ++i) {
        switch (s[i]) {
            case '&':
                out += "&amp;";
                break;
            case '<':
                out += "&lt;";
                break;
            case '>':
                out += "&gt;";
                break;
            case '"':
                out += "&quot;";
                break;
            default:
                out += s[i];
                break;
        }
    }
    return out;
}

void write_csv(const std::string& path, const std::vector<Series>& series) {
    std::ofstream out(path.c_str());
    out << std::setprecision(15);
    out << "series,x,y\n";
    for (std::size_t s = 0; s < series.size(); ++s) {
        for (std::size_t i = 0; i < series[s].values.size(); ++i) {
            out << series[s].name << "," << series[s].values[i].x << ","
                << series[s].values[i].y << "\n";
        }
    }
}

void write_svg(const std::string& path,
               const std::string& title,
               const std::string& xlabel,
               const std::string& ylabel,
               const std::vector<Series>& series) {
    const int width = 860;
    const int height = 560;
    const int left = 82;
    const int right = 30;
    const int top = 48;
    const int bottom = 72;

    double xmin = 0.0;
    double xmax = 1.0;
    double ymin = 0.0;
    double ymax = 1.0;
    bool first = true;

    for (std::size_t s = 0; s < series.size(); ++s) {
        for (std::size_t i = 0; i < series[s].values.size(); ++i) {
            const double x = series[s].values[i].x;
            const double y = series[s].values[i].y;
            if (!std::isfinite(x) || !std::isfinite(y)) {
                continue;
            }
            if (first) {
                xmin = xmax = x;
                ymin = ymax = y;
                first = false;
            } else {
                xmin = std::min(xmin, x);
                xmax = std::max(xmax, x);
                ymin = std::min(ymin, y);
                ymax = std::max(ymax, y);
            }
        }
    }

    if (std::fabs(xmax - xmin) < 1.0e-14) {
        xmax = xmin + 1.0;
    }
    if (std::fabs(ymax - ymin) < 1.0e-14) {
        ymax = ymin + 1.0;
    }

    const double xpad = 0.04 * (xmax - xmin);
    const double ypad = 0.08 * (ymax - ymin);
    xmin -= xpad;
    xmax += xpad;
    ymin -= ypad;
    ymax += ypad;
    if (ymin > 0.0) {
        ymin = 0.0;
    }

    const double plot_w = static_cast<double>(width - left - right);
    const double plot_h = static_cast<double>(height - top - bottom);

    auto sx = [&](double x) {
        return left + (x - xmin) / (xmax - xmin) * plot_w;
    };
    auto sy = [&](double y) {
        return top + (ymax - y) / (ymax - ymin) * plot_h;
    };

    std::ofstream out(path.c_str());
    out << std::fixed << std::setprecision(3);
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width
        << "\" height=\"" << height << "\" viewBox=\"0 0 " << width << " "
        << height << "\">\n";
    out << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
    out << "<text x=\"" << width / 2
        << "\" y=\"28\" text-anchor=\"middle\" font-family=\"Arial\""
        << " font-size=\"18\" font-weight=\"700\">" << escape_xml(title)
        << "</text>\n";
    out << "<line x1=\"" << left << "\" y1=\"" << height - bottom
        << "\" x2=\"" << width - right << "\" y2=\"" << height - bottom
        << "\" stroke=\"#222\" stroke-width=\"1.5\"/>\n";
    out << "<line x1=\"" << left << "\" y1=\"" << top << "\" x2=\"" << left
        << "\" y2=\"" << height - bottom
        << "\" stroke=\"#222\" stroke-width=\"1.5\"/>\n";

    for (int tick = 0; tick <= 5; ++tick) {
        const double x = xmin + (xmax - xmin) * tick / 5.0;
        const double px = sx(x);
        out << "<line x1=\"" << px << "\" y1=\"" << height - bottom
            << "\" x2=\"" << px << "\" y2=\"" << height - bottom + 6
            << "\" stroke=\"#222\"/>\n";
        out << "<text x=\"" << px << "\" y=\"" << height - bottom + 24
            << "\" text-anchor=\"middle\" font-family=\"Arial\""
            << " font-size=\"12\">" << std::setprecision(4) << x
            << std::setprecision(3) << "</text>\n";

        const double y = ymin + (ymax - ymin) * tick / 5.0;
        const double py = sy(y);
        out << "<line x1=\"" << left - 6 << "\" y1=\"" << py
            << "\" x2=\"" << left << "\" y2=\"" << py
            << "\" stroke=\"#222\"/>\n";
        out << "<line x1=\"" << left << "\" y1=\"" << py << "\" x2=\""
            << width - right << "\" y2=\"" << py
            << "\" stroke=\"#e8e8e8\"/>\n";
        out << "<text x=\"" << left - 10 << "\" y=\"" << py + 4
            << "\" text-anchor=\"end\" font-family=\"Arial\""
            << " font-size=\"12\">" << std::setprecision(4) << y
            << std::setprecision(3) << "</text>\n";
    }

    out << "<text x=\"" << width / 2 << "\" y=\"" << height - 18
        << "\" text-anchor=\"middle\" font-family=\"Arial\""
        << " font-size=\"15\">" << escape_xml(xlabel) << "</text>\n";
    out << "<text transform=\"translate(24 " << height / 2
        << ") rotate(-90)\" text-anchor=\"middle\" font-family=\"Arial\""
        << " font-size=\"15\">" << escape_xml(ylabel) << "</text>\n";

    for (std::size_t s = 0; s < series.size(); ++s) {
        if (series[s].values.empty()) {
            continue;
        }
        out << "<polyline fill=\"none\" stroke=\"" << series[s].color
            << "\" stroke-width=\"2.2\" points=\"";
        for (std::size_t i = 0; i < series[s].values.size(); ++i) {
            const double x = series[s].values[i].x;
            const double y = series[s].values[i].y;
            if (std::isfinite(x) && std::isfinite(y)) {
                out << sx(x) << "," << sy(y) << " ";
            }
        }
        out << "\"/>\n";
        if (series[s].points) {
            for (std::size_t i = 0; i < series[s].values.size(); ++i) {
                const double x = series[s].values[i].x;
                const double y = series[s].values[i].y;
                if (std::isfinite(x) && std::isfinite(y)) {
                    out << "<circle cx=\"" << sx(x) << "\" cy=\"" << sy(y)
                        << "\" r=\"2.6\" fill=\"" << series[s].color
                        << "\"/>\n";
                }
            }
        }
    }

    int legend_y = top + 14;
    for (std::size_t s = 0; s < series.size(); ++s) {
        const int lx = width - right - 185;
        out << "<line x1=\"" << lx << "\" y1=\"" << legend_y
            << "\" x2=\"" << lx + 28 << "\" y2=\"" << legend_y
            << "\" stroke=\"" << series[s].color
            << "\" stroke-width=\"2.5\"/>\n";
        out << "<text x=\"" << lx + 36 << "\" y=\"" << legend_y + 4
            << "\" font-family=\"Arial\" font-size=\"12\">"
            << escape_xml(series[s].name) << "</text>\n";
        legend_y += 18;
    }

    out << "</svg>\n";
}

}  // namespace

void ensure_dir(const std::string& path) {
    if (::mkdir(path.c_str(), 0775) != 0 && errno != EEXIST) {
        throw std::runtime_error("could not create directory: " + path);
    }
}

Series make_series(const std::string& name,
                   const std::string& color,
                   bool points,
                   const std::vector<Point>& values) {
    Series s;
    s.name = name;
    s.color = color;
    s.points = points;
    s.values = values;
    return s;
}

std::vector<Point> points_from_arrival(const ArrivalResult& r,
                                       const std::vector<double>& y) {
    std::vector<Point> pts;
    for (int i = 1; i <= r.n; ++i) {
        Point p;
        p.x = r.t[i];
        p.y = y[i];
        pts.push_back(p);
    }
    return pts;
}

void write_plot(const std::string& out_dir,
                const std::string& stem,
                const std::string& title,
                const std::string& xlabel,
                const std::string& ylabel,
                const std::vector<Series>& series) {
    write_csv(out_dir + "/" + stem + ".csv", series);
    write_svg(out_dir + "/" + stem + ".svg", title, xlabel, ylabel, series);
}

}  // namespace qtoa
