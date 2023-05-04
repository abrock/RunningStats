#undef NDEBUG
#include <cassert>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

#include "runningstats.h"

#include "gnuplot-iostream.h"

#include "misc.h"

namespace runningstats {

Histogram2Dfixed::Histogram2Dfixed(
        const double _width_1,
        const double _width_2,
        const double _min_1,
        const double _min_2,
        const double _max_1,
        const double _max_2) :
    width_1(sanitize_bin_width(_width_1, _min_1, _max_1)),
    width_2(sanitize_bin_width(_width_2, _min_2, _max_2)),
    min_1(_min_1),
    min_2(_min_2),
    max_1(_max_1),
    max_2(_max_2){
    double const range_1 = std::abs(max_1 - min_1);
    double const range_2 = std::abs(max_2 - min_2);
    data = std::vector<std::vector<size_t>>
                                            (std::ceil(range_1 / width_1)+1,
                                             std::vector<size_t>(std::ceil(range_2 / width_2)+1, 0));
    stats_per_column = std::vector<QuantileStats<float> >(std::ceil(range_1 / width_1)+1);
}

Histogram2Dfixed::Histogram2Dfixed(const Histogram2Dfixed &rhs):
    width_1(rhs.width_1),
    width_2(rhs.width_2),
    min_1(rhs.min_1),
    min_2(rhs.min_2),
    max_1(rhs.max_1),
    max_2(rhs.max_2),
    data(rhs.data),
    stats_per_column(rhs.stats_per_column),
    total_count(rhs.total_count)
{}

Histogram2Dfixed::~Histogram2Dfixed() {
    data.clear();
}

double Histogram2Dfixed::sanitize_bin_width(const double width, const double min, const double max) {
    if (!std::isfinite(width) || !std::isfinite(min) || !std::isfinite(max)) {
        return 1;
    }
    if (width <= 0
            || (std::abs(max) + std::abs(min)) / width >= std::vector<double>().max_size()) {
        return std::max(max-min, 1.0);
    }
    return width;
}

bool Histogram2Dfixed::push_unsafe(const double val1, const double val2) {
    if (!std::isfinite(val1) || !std::isfinite(val2) || val1 < min_1 || val1 > max_1 || val2 < min_2 || val2 > max_2) {
        return false;
    }
    size_t const col = std::round((val1 - min_1) / width_1);
    size_t const row = std::round((val2 - min_2) / width_2);
    data[col][row]++;
    stats_per_column[col].push_unsafe(val2);
    total_count++;
    return true;
}

void Histogram2Dfixed::plotHist(const std::string prefix, HistConfig const& conf) const {
    std::string const data_file = prefix + ".data";
    Misc::make_target_dir(data_file);
    std::ofstream data_out(data_file);
    if (conf.normalize_x) {
        for (size_t row = 0; row < data.size(); ++row) {
            std::vector<size_t> const& row_data = data[row];
            size_t row_sum = 0;
            for (size_t col = 0; col < row_data.size(); ++col) {
                row_sum += row_data[col];
            }
            if (row_sum < 20) {
                row_sum = 20;
            }
            double const row_bin = width_1 * row + min_1;
            for (size_t col = 0; col < row_data.size(); ++col) {
                double const col_bin = width_2 * col + min_2;
                data_out << row_bin << "\t" << col_bin << "\t" << (double(row_data[col]) / (row_sum * width_1 * width_2)) << std::endl;
            }
            data_out << std::endl;
        }
    }
    else {
        for (size_t row = 0; row < data.size(); ++row) {
            std::vector<size_t> const& row_data = data[row];
            double const row_bin = width_1 * row + min_1;
            for (size_t col = 0; col < row_data.size(); ++col) {
                double const col_bin = width_2 * col + min_2;
                data_out << row_bin << "\t" << col_bin << "\t" << (conf.absolute ? row_data[col] : double(row_data[col]) / (total_count * width_1 * width_2)) << std::endl;
            }
            data_out << std::endl;
        }
    }
    std::stringstream cmd;
    bool const range_1_empty = (min_1 == max_1);
    bool const range_2_empty = (min_2 == max_2);
    cmd << "#!/usr/bin/gnuplot \n";
    cmd << "set term png;\n";
    cmd << "set output '" << prefix << ".png';\n";
    cmd << conf.toString();
    cmd << "set xrange[" << min_1 - (range_1_empty ? 1:width_1/2) << " : " << max_1 + (range_1_empty ? 1:width_1/2) << "];\n";
    cmd << "set yrange[" << min_2 - (range_2_empty ? 1:width_2/2) << " : " << max_2 + (range_2_empty ? 1:width_2/2) << "];\n";
    cmd << "set xtics out;\n";
    cmd << "set ytics out;\n";
    cmd << conf.colorMapCmd();
    cmd << "plot '" << data_file << "' u 1:2:3 with image notitle ";
    cmd << conf.plotContours(data_file);
    cmd << conf.plotLines(data_file);
    cmd << ";\n";
    cmd << conf.generateFormatCommands(prefix);
    //cmd << "set term tikz;\n";
    //cmd << "set output '" << prefix << ".tex';\n";
    //cmd << "replot;\n";
    std::ofstream cmd_out(prefix + ".gpl");
    cmd_out << cmd.str();
    cmd_out.close();
    gnuplotio::Gnuplot plt;
    plt << cmd.str();
}

void Histogram2Dfixed::plotHistPm3D(const std::string prefix, const HistConfig &conf) const {
    std::string const data_file = prefix + ".data";
    Misc::make_target_dir(data_file);
    std::ofstream data_out(data_file);
    if (conf.normalize_x) {
        size_t col_sum = 0;
        for (size_t col = 0; col < data.size(); ++col) {
            std::vector<size_t> const& col_data = data[col];
            col_sum = 0;
            for (size_t row = 0; row < col_data.size(); ++row) {
                col_sum += col_data[row];
            }
            if (col_sum < 20) {
                col_sum = 20;
            }
            double const row_bin = width_1 * col + min_1 - width_1/2;
            for (size_t row = 0; row < col_data.size(); ++row) {
                double const col_bin = width_2 * row + min_2 - width_2/2;
                data_out << row_bin << "\t" << col_bin << "\t" << (double(col_data[row]) / (col_sum * width_1 * width_2)) << std::endl;
            }
            data_out << row_bin << "\t" << width_2 * col_data.size() + min_2 - width_2/2
                     << "\t" << (double(col_data.back()) / (col_sum * width_1 * width_2)) << std::endl;
            data_out << std::endl;
        }
        auto const& col_data = data.back();
        double const col_bin = width_1 * data.size() + min_1 - width_1/2;
        for (size_t row = 0; row < col_data.size(); ++row) {
            double const row_bin = width_2 * row + min_2 - width_2/2;
            data_out << col_bin << "\t" << row_bin << "\t" << (double(col_data[row]) / (col_sum * width_1 * width_2)) << std::endl;
        }
        data_out << col_bin << "\t" << width_2 * col_data.size() + min_2 - width_2/2
                 << "\t" << (double(col_data.back()) / (col_sum * width_1 * width_2)) << std::endl;
        data_out << std::endl;
    }
    else {
        std::string row_text;
        for (size_t row = 0; row < data.size(); ++row) {
            std::vector<size_t> const& row_data = data[row];
            double const row_bin = width_1 * row + min_1 - width_1/2;
            for (size_t col = 0; col < row_data.size(); ++col) {
                double const col_bin = width_2 * col + min_2 - width_2/2;
                data_out << row_bin << "\t" << col_bin << "\t" << (conf.absolute ? row_data[col] : double(row_data[col]) / (total_count * width_1 * width_2)) << std::endl;
            }
            data_out << row_bin << "\t" << width_2 * row_data.size() + min_2 - width_2/2 << "\t" << (conf.absolute ? row_data.back() : double(row_data.back()) / (total_count * width_1 * width_2)) << std::endl;
            data_out << std::endl;
        }
        std::vector<size_t> const& row_data = data.back();
        double const row_bin = width_1 * data.size() + min_1 - width_1/2;
        for (size_t col = 0; col < row_data.size(); ++col) {
            double const col_bin = width_2 * col + min_2 - width_2/2;
            data_out << row_bin << "\t" << col_bin << "\t" << (conf.absolute ? row_data[col] : double(row_data[col]) / (total_count * width_1 * width_2)) << std::endl;
        }
        data_out << row_bin << "\t" << width_2 * row_data.size() + min_2 - width_2/2 << "\t" << (conf.absolute ? row_data.back() : double(row_data.back()) / (total_count * width_1 * width_2)) << std::endl;
        data_out << std::endl;
    }
    std::stringstream cmd;
    gnuplotio::Gnuplot plt;
    bool const range_1_empty = (min_1 == max_1);
    bool const range_2_empty = (min_2 == max_2);
    cmd << "#!/usr/bin/gnuplot \n";
    cmd << "set term png;\n";
    cmd << "set output '" << prefix << ".png';\n";
    cmd << conf.toString();
    cmd << "set xrange[" << min_1 - (range_1_empty ? 1:width_1/2) << " : " << max_1 + (range_1_empty ? 1:width_1/2) << "];\n";
    cmd << "set yrange[" << min_2 - (range_2_empty ? 1:width_2/2) << " : " << max_2 + (range_2_empty ? 1:width_2/2) << "];\n";
    cmd << "set xtics out;\n";
    cmd << "set ytics out;\n";
    cmd << "set view map;\n";
    cmd << "set pm3d corners2color c1;\n";
    cmd << "set key out horiz;\n";
    cmd << conf.colorMapCmd();
    cmd << "splot '" << data_file << "' u 1:2:3 with pm3d notitle";
    for (HistConfig::ExtractConf const& extractor : conf.extractors) {
        std::vector<std::pair<double, double> > data_out;
        data_out.reserve(data.size());
        for (size_t col = 0; col < data.size(); ++col) {
            QuantileStats<float> const& col_data = stats_per_column[col];
            double const row_bin = width_1 * col + min_1;
            data_out.push_back({row_bin, col_data.getStat(extractor.ex, extractor.value)});
        }
        std::string extract_name = extractor.getName();
        std::string title = extractor.title.empty() ? " notitle " : " title '" + escape(extractor.title) + "' ";
        std::string color = extractor.color.empty() ? " " : " lc rgb '" + escape(extractor.color) + "' ";
        cmd << ", " << plt.file(data_out, prefix + extract_name + ".data") << " u 1:2:(0.0) w l " << color << title;
    }
    cmd << conf.plotContours(data_file);
    cmd << conf.plotLines(data_file);
    cmd << ";\n";
    cmd << conf.generateFormatCommands(prefix);
    //cmd << "set term png;\n";
    //cmd << "set output '" << prefix << ".png';\n";
    //cmd << "replot;\n";
    //cmd << "set term tikz;\n";
    //cmd << "set output '" << prefix << ".tex';\n";
    //cmd << "replot;\n";
    std::ofstream cmd_out(prefix + ".gpl");
    cmd_out << cmd.str();
    cmd_out.close();
    plt << cmd.str();
}

void Histogram2Dfixed::plotHist(const std::string prefix, const bool absolute) const {
    HistConfig conf;
    conf.absolute = absolute;
    plotHist(prefix, conf);
}

} // namespace runningstats
