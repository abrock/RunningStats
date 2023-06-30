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

Histogram::Histogram(const double _bin_size) : bin_size(_bin_size) {

}

Histogram::Histogram(const Histogram &rhs): bin_size(rhs.bin_size), data(rhs.data) {

}

bool Histogram::push(double value) {
    if (!std::isfinite(value)) {
        return false;
    }
    std::lock_guard guard(push_mutex);
    push_unsafe(value);
    return true;
}

template<class T>
bool Histogram::push_vector(const std::vector<T> &values) {
    std::lock_guard guard(push_mutex);
    return push_vector_unsafe(values);
}

bool Histogram::push_unsafe(const double value) {
    if (!std::isfinite(value)) {
        return false;
    }
    int64_t const bin = int64_t(std::round(value/bin_size));
    data[bin]++;
    RunningStats::push_unsafe(value);

    return true;
}

template<class T>
bool Histogram::push_vector_unsafe(const std::vector<T> &values) {
    bool success = true;
    for (const double  val : values) {
        success &= push_unsafe(val);
    }
    return success;
}

template bool Histogram::push_vector_unsafe(const std::vector<unsigned long>&);
template bool Histogram::push_vector_unsafe(const std::vector<float>&);
template bool Histogram::push_vector_unsafe(const std::vector<double>&);


std::vector<std::pair<double, double> > Histogram::getAbsoluteHist() const {
    std::vector<std::pair<double, double> > result;
    if (data.empty()) {
        return result;
    }
    for (int64_t ii = data.begin()->first; ii <= data.rbegin()->first; ++ii) {
        auto const it = data.find(ii);
        if (it != data.end()) {
            result.push_back(std::pair<double, double>(it->first * bin_size, double(it->second)));
        }
        else {
            result.push_back(std::pair<double, double>(ii * bin_size, 0.0));
        }
    }
    return result;
}

std::vector<std::pair<double, double> > Histogram::getRelativeHist() const {
    std::vector<std::pair<double, double> > result;
    if (data.empty()) {
        return result;
    }
    for (int64_t ii = data.begin()->first; ii <= data.rbegin()->first; ++ii) {
        auto const it = data.find(ii);
        if (it != data.end()) {
            result.push_back(std::pair<double, double>(it->first * bin_size, double(it->second)/(double(n)*bin_size)));
        }
        else {
            result.push_back(std::pair<double, double>(ii * bin_size, 0.0));
        }
    }
    return result;
}

void Histogram::plotHist(const std::string prefix, const HistConfig conf) const {
    if (data.empty()) {
        return;
    }
    gnuplotio::Gnuplot gpl;
    std::stringstream cmd;
    disable_thousands_separator(cmd);
    std::string data_file = prefix + ".data";
    Misc::make_target_dir(data_file);

    cmd << "#!/usr/bin/gnuplot \n";
    cmd << "set term png;\n";
    cmd << "set output '" << prefix + ".png'; \n";
    cmd << conf.toString();
    cmd << "set title '" << escape(conf.title) << " n=" << escape(tostring(getCount())) << ", m=" << getMean() << ", s=" << getStddev() << "';\n";
    if (!conf.dataLabel.empty()) {
        cmd << "set xlabel '" << escape(conf.dataLabel) << "';\n";
    }
    cmd << "set ylabel 'Estimated PDF';\n";

    cmd << getXrange(conf);

    auto const hist_data = conf.absolute ? getAbsoluteHist() : getRelativeHist();

    cmd << "plot " << gpl.file1d(hist_data, data_file) << " w boxes notitle; \n";
    //cmd << "set term tikz; \n";
    //cmd << "set output '" << prefix << ".tex'; \n";
    //cmd << "replot;\n";
    cmd << conf.generateFormatCommands(prefix);

    gpl << cmd.str();

    std::ofstream out(prefix + ".gpl");
    out << cmd.str();
}

std::string Histogram::printHist() const {
    std::stringstream out;
    std::vector<std::pair<double, double> > hist = getAbsoluteHist();
    for (size_t ii = 0; ii < hist.size(); ++ii) {
        out << ii << ": " << hist[ii].first << ":\t" << hist[ii].second << "\t(" << 100.0*hist[ii].second/n << "%)" << std::endl;
    }
    return out.str();
}

void Histogram::plotHist(const std::string prefix, const bool absolute) const {
    plotHist(prefix, HistConfig().setAbsolute(absolute));
}

size_t Histogram::getBinCount(const double value) const {
    int64_t const bin = int64_t(std::round(value/bin_size));
    auto const it = data.find(bin);
    if (it != data.end()) {
        return it->second;
    }
    return 0;
}

double Histogram::getLikelihood() const {
    double result = 1;
    for (std::pair<const int64_t, size_t> const& it : data) {
        double const pi_k = double(it.second) / double(n);
        double const p_k = pi_k / bin_size;
        result *= p_k;
    }
    return result;
}

double Histogram::getPosteriorProbability() const {
    double result = 0;
    double const N = n; // Total number of samples
    double const M = (getMax() - getMin()) / bin_size; // Total number of bins (including empty ones)
    result += n * std::log(M)
            + std::lgamma(M/2)
            - M * std::lgamma(0.5)
            - std::lgamma(N + M/2);
    int64_t const start = data.begin()->first;
    int64_t const end = data.rbegin()->first;
    for (int64_t ii = start; ii <= end; ++ii) {
        auto const it = data.find(ii);
        if (it != data.end()) {
            result += std::lgamma(.5 + double(it->second));
        }
        else {
            result += std::lgamma(.5);
        }
    }
    return result;
}

std::string Histogram::getXrange(const HistConfig &conf) const {
    double _min = std::max(min, conf.min_x);
    double _max = std::min(max, conf.max_x);
    std::stringstream out;
    disable_thousands_separator(out);
    out << "set xrange[" << _min << ":" << _max << "];" << std::endl;
    return out.str();
}

Histogram2D::Histogram2D(const double _bin_1, const double _bin_2) : width_1(_bin_1), width_2(_bin_2) {
}

Histogram2D::Histogram2D(const Histogram2D &rhs) : width_1(rhs.width_1), width_2(rhs.width_2), data(rhs.data) {}

bool Histogram2D::push(const double val1, const double val2) {
    if (!std::isfinite(val1) || !std::isfinite(val2)) {
        return false;
    }
    const std::lock_guard<std::mutex> lock(push_mutex);
    push_unsafe(val1, val2);
    return true;
}

bool Histogram2D::push_unsafe(const double val1, const double val2) {
    if (!std::isfinite(val1) || !std::isfinite(val2)) {
        return false;
    }
    int64_t const bin1 = int64_t(std::round(val1/width_1));
    int64_t const bin2 = int64_t(std::round(val2/width_2));
    data[bin1][bin2]++;
    bins_1.insert(bin1);
    bins_2.insert(bin2);
    total_count++;
    return true;
}

void Histogram2D::plotHist(const std::string prefix, const bool absolute) const {
    std::string const data_file = prefix + ".data";
    Misc::make_target_dir(data_file);
    std::ofstream data_out(data_file);
    disable_thousands_separator(data_out);
    for (auto const& it : data) {
        double const row = width_1 * double(it.first);
        for (const double col_id: bins_2) {
            double const col = width_2 * double(col_id);
            auto const found = it.second.find(col_id);
            if (it.second.end() != found) {
                size_t const count = found->second;
                double const probability = double(count)/double(total_count);
                double const density = probability / (width_1 * width_2);
                double const val = absolute ? count : density;
                data_out << row << "\t" << col << "\t" << val << std::endl;
            }
            else {
                data_out << row << "\t" << col << "\t" << 0 <<  std::endl;
            }
        }
        /*
        for (auto const& it2 : it.second) {
            double const col = width_2 * double(it2.first);
            size_t const count = it2.second;
            double const probability = double(count)/double(total_count);
            double const density = probability / (width_1 * width_2);
            double const val = absolute ? count : density;
            data_out << row << "\t" << col << "\t" << val << std::endl;
        }
        */
        data_out << std::endl;
    }
    std::stringstream cmd;
    disable_thousands_separator(cmd);
    cmd << "#!/usr/bin/gnuplot \n";
    cmd << "set term png;\n";
    cmd << "set output '" << prefix << ".png';\n";
    cmd << "plot '" << data_file << "' u 2:1:3 with image notitle;\n";
    gnuplotio::Gnuplot plt;
    plt << cmd.str();
    std::ofstream cmd_out(prefix + ".gpl");
    cmd_out << cmd.str();
}

} // namespace runningstats
