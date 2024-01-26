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

template<class T>
T QuantileStats<T>::getMedian() const {
    return getQuantile(0.5);
}

template<class T>
void QuantileStats<T>::push(const double value){
    std::lock_guard guard(push_mutex);
    push_unsafe(value);
}

template<class T>
void QuantileStats<T>::push_unsafe(const double value) {
    sorted = false;
    values.push_back(value);
    RunningStats::push_unsafe(value);
}

template<class T>
void QuantileStats<T>::push_repeated(const double value, const size_t count) {
    if (!std::isfinite(value)) {
        return;
    }
    std::lock_guard guard(push_mutex);
    push_repeated_unsafe(value, count);
}

template<class T>
void QuantileStats<T>::push_repeated_unsafe(const double value, const size_t count) {
    if (!std::isfinite(value)) {
        return;
    }
    for (size_t ii = 0; ii < count; ++ii) {
        push_unsafe(value);
    }
}

template<class T>
template<class U>
void QuantileStats<T>::push_unsafe(QuantileStats<U> const& other) {
    push_unsafe(other.values);
}

template void QuantileStats<float>::push_unsafe(QuantileStats<float> const& other);

template<class T>
Histogram QuantileStats<T>::getHistogram(const double bin_size) {
    Histogram result(bin_size);
    result.push_vector_unsafe(values);
    return result;
}


template<class T>
T QuantileStats<T>::getQuantile(const double quantile) const {
    if (quantile <= 0) {
        return min;
    }
    if (quantile >= 1) {
        return max;
    }
    return getQuantile(quantile, values, sorted);
}

template<class T>
T QuantileStats<T>::getQuantile(const double quantile, std::vector<T>& values, bool & sorted)  {
    if (quantile <= 0) {
        return getMin(values);
    }
    if (quantile >= 1) {
        return getMax(values);
    }
    if (values.size() == 0) {
        return 0;
    }
    if (values.size() == 1) {
        return values[0];
    }
    //sort();
    size_t const n = static_cast<size_t>(quantile * (values.size()-1));
    if (sorted) {
        return values[n];
    }
    std::nth_element(values.begin(), values.begin() + n, values.end());
    sorted = false;
    return values[n];
}

template<class T>
double QuantileStats<T>::getMin(const std::vector<T> &values) {
    if (values.empty()) {
        return 0;
    }
    T result = values[0];
    for (size_t ii = 1; ii < values.size(); ++ii) {
        if (values[ii] < result) {
            result = values[ii];
        }
    }
    return double(result);
}

template<class T>
double QuantileStats<T>::getMax(const std::vector<T> &values) {
    if (values.empty()) {
        return 0;
    }
    T result = values[0];
    for (size_t ii = 1; ii < values.size(); ++ii) {
        if (values[ii] > result) {
            result = values[ii];
        }
    }
    return double(result);
}

template<class T>
std::string QuantileStats<T>::getXrange(const HistConfig &conf) const {
    double _min = std::max(min, conf.min_x);
    double _max = std::min(max, conf.max_x);

    if (conf.ignore_amount > 0) {
        _min = std::max(_min, double(getQuantile(conf.ignore_amount/2)));
        _max = std::min(_max, double(getQuantile(1.0 - conf.ignore_amount/2)));
    }
    std::stringstream out;
    disable_thousands_separator(out);
    out << "set xrange[" << _min << ":" << _max << "];\n";
    return out.str();
}

template<class T>
void QuantileStats<T>::setRangeByIgnoreAmount(HistConfig &conf) const {
    if (conf.ignore_amount > 0) {
        conf.min_x = std::max(conf.min_x, double(getQuantile(conf.ignore_amount/2)));
        conf.max_x = std::min(conf.max_x, double(getQuantile(1.0 - conf.ignore_amount/2)));
    }
}

template<class T>
void QuantileStats<T>::saveSummary(const std::string &filename) {
    std::ofstream out(filename);
    getSummary(out);
}

template<class T>
std::string QuantileStats<T>::getSummary() {
    std::stringstream out;
    getSummary(out);
    return out.str();
}

template<class T>
void QuantileStats<T>::getSummary(std::ostream &out) {
    out << "Quantiles 1%, 5%, 10%, 15%, 25%, 50%, 75%, 85%, 90%, 95%, 99%:\n";
    for (double q : {0.01, 0.05, 0.1, 0.15, 0.25, 0.5, 0.75, 0.85, 0.9, 0.95, 0.99}) {
        out << getQuantile(q) << "\t";
    }
    out << "\n" << "Median, IQR:\n" << getQuantile(0.5) << " [ " << getQuantile(0.25) << " " << getQuantile(0.75) << "\n"
        << "Mean, stddev: " << getMean() << "\t" << getStddev() << std::endl;
}

template<class T>
double QuantileStats<T>::getStat(const HistConfig::Extract type, const double param) const {
    switch (type) {
    case HistConfig::Extract::Mean: return getMean();
    case HistConfig::Extract::Median: return getMedian();
    case HistConfig::Extract::Quantile: return getQuantile(param);
    case HistConfig::Extract::Stddev: return getStddev();
    case HistConfig::Extract::TrimmedMean: return getTrimmedMean(param);
    case HistConfig::Extract::Variance: return getVar();
    case HistConfig::Extract::Count: return getCount();
    default: break;
    }
    throw std::runtime_error("Extract type not available");
}

template<class T>
double QuantileStats<T>::getStat(const std::pair<HistConfig::Extract, double> pair) const {
    return getStat(pair.first, pair.second);
}

template<class T>
void QuantileStats<T>::save(std::ostream &out) const {
    std::lock_guard guard(push_mutex);
    for (size_t ii = 0; ii < values.size(); ++ii) {
        WRITE_BIN(out, values[ii]);
    }
}

template<class T>
void QuantileStats<T>::save(const std::string &filename) const {
    std::ofstream out(filename);
    save(out);
}

template<class T>
void QuantileStats<T>::load(std::istream &in) {
    std::lock_guard guard(push_mutex);
    sorted = false;
    while (in) {
        T value = 0;
        if (READ_BIN(in, value)) {
            push_unsafe(value);
        }
    }
}

template<class T>
void QuantileStats<T>::load(const std::string &filename) {
    std::ifstream in(filename);
    load(in);
}

template<class T>
T QuantileStats<T>::getInverseQuantile(const double value) const {
    if (value <= min) {
        return 0;
    }
    if (value >= max) {
        return 1;
    }
    sort();

    typename std::vector<T>::iterator low = std::lower_bound (values.begin(), values.end(), value); //          ^
    typename std::vector<T>::iterator up = std::upper_bound (values.begin(), values.end(), value); //          ^

    return (static_cast<double>(low - values.begin()) + static_cast<double>(up - values.begin())) / (2*values.size());
}

template<class T>
void QuantileStats<T>::sort() const {
    if (!sorted) {
        std::lock_guard guard(push_mutex);
        std::sort(values.begin(), values.end());
        sorted = true;
    }
}

template<class T>
void QuantileStats<T>::plotHist(const std::string prefix, const double bin_size, HistConfig conf) const {
    double _bin_size = bin_size > 0 ? bin_size : FreedmanDiaconisBinSize();
    setRangeByIgnoreAmount(conf);
    double x_range = getMax() - getMin();
    if (conf.ignore_amount > 0) {
        x_range = conf.max_x - conf.min_x;
    }
    double const n_x = x_range / _bin_size;
    if (std::isfinite(n_x) && conf.max_nx > 0 && conf.max_nx < n_x) {
        _bin_size = x_range / conf.max_nx;
    }
    if (x_range / _bin_size > conf.max_plot_pts && conf.max_plot_pts > 0) {
        _bin_size = x_range / conf.max_plot_pts;
    }
    Histogram h(_bin_size);
    h.push_vector_unsafe(values);
    h.plotHist(prefix, conf);
}

template<class T>
void QuantileStats<T>::plotHist(const std::string prefix, HistConfig conf) const {
    plotHist(prefix, FreedmanDiaconisBinSize(), conf);
}

template<class T>
void QuantileStats<T>::plotHist(const std::string prefix, const double bin_size, const bool absolute) const {
    double const _bin_size = bin_size > 0 ? bin_size : FreedmanDiaconisBinSize();
    Histogram h(_bin_size);
    h.push_vector_unsafe(values);
    h.plotHist(prefix, HistConfig().setAbsolute(absolute));
}

namespace  {
size_t const max_pts_CDF_plot = 5000;
}

template<class T>
void QuantileStats<T>::plotCDF(const std::string prefix, HistConfig conf) const {
    if (values.size() < 2) {
        return;
    }
    if (values.size() > max_pts_CDF_plot) {
        plotReducedCDF(prefix, conf);
        return;
    }
    sort();

    gnuplotio::Gnuplot gpl;
    std::stringstream cmd;
    disable_thousands_separator(cmd);
    std::string data_file = prefix + ".data";
    Misc::make_target_dir(data_file);
    cmd << "#!/usr/bin/gnuplot \n";
    cmd << "set term png;\n"
        << "set output '" << prefix + ".png'; \n"
        << conf.toString()
        << "set title '" << escape(conf.title) << " n=" << escape(tostring(getCount())) << ", m=" << getMean() << ", s=" << getStddev() << "';\n";

    if (!conf.dataLabel.empty()) {
        cmd << "set xlabel '" << escape(conf.dataLabel) << "';\n";
    }
    cmd << "set ylabel 'Estimated CDF';\n";


    cmd << getXrange(conf);

    {
        std::lock_guard guard(push_mutex);
        cmd << "plot " << gpl.file(values, data_file) << " u 1:"
            << (conf.absolute ? "0" : "($0/" + std::to_string(values.size()-1) + ")") << " w l notitle; \n";
    }
    //cmd << "set term tikz; \n"
    //<< "set output '" << prefix << ".tex'; \n"
    //<< "replot;\n";
    cmd << conf.generateFormatCommands(prefix);

    gpl << cmd.str();

    std::ofstream out(prefix + ".gpl");
    out << cmd.str();
}

template<class T>
void QuantileStats<T>::plotReducedCDF(const std::string prefix, HistConfig conf) const {
    if (values.size() < 2) {
        return;
    }
    if (values.size() < max_pts_CDF_plot) {
        plotCDF(prefix, conf);
        return;
    }
    sort();
    T current = min;
    std::vector<std::pair<T, float> > plot_values;
    plot_values.reserve(max_pts_CDF_plot);
    double const _min = std::max(min, conf.min_x);
    double const _max = std::min(max, conf.max_x);
    double const orig_range = max - min;
    double const plot_range = _max - _min;
    double reduced_range_factor = orig_range / plot_range;
    if (!std::isfinite(reduced_range_factor) || reduced_range_factor < 1) {
        reduced_range_factor = 1;
    }
    {
        std::lock_guard guard(push_mutex);
        for (size_t ii = 0; ii < values.size(); ++ii) {
            if (values[ii] >= current) {
                plot_values.push_back({values[ii], float(ii)/(conf.absolute ? 1 : values.size()-1)});
                current = min + ((max-min) * plot_values.size()) / (max_pts_CDF_plot * reduced_range_factor);
            }
        }
    }
    plot_values.push_back({max, conf.absolute ? values.size() : 1});
    {
        auto prev = plot_values.front();
        for (size_t ii = 0; ii < plot_values.size(); ++ii) {
            auto const& it = plot_values[ii];
            assert(it.first >= prev.first);
            assert(it.second >= prev.second);
            prev = it;
        }
    }

    gnuplotio::Gnuplot gpl;
    std::stringstream cmd;
    disable_thousands_separator(cmd);
    std::string data_file = prefix + ".data";
    Misc::make_target_dir(data_file);
    cmd << "#!/usr/bin/gnuplot \n";
    cmd << "set term png;\n"
        << "set output '" << prefix + ".png'; \n"
        << conf.toString()
        << "set title '" << escape(conf.title) << " n=" << escape(tostring(getCount())) << ", m=" << getMean() << ", s=" << getStddev() << "';\n";

    if (!conf.dataLabel.empty()) {
        cmd << "set xlabel '" << escape(conf.dataLabel) << "';\n";
    }
    cmd << "set ylabel 'Estimated CDF';\n";

    cmd << getXrange(conf);

    cmd << "plot " << gpl.file(plot_values, data_file) << " u 1:2 w l notitle; \n";
    //cmd << "set term tikz; \n"
    // << "set output '" << prefix << ".tex'; \n"
    // << "replot;\n";
    cmd << conf.generateFormatCommands(prefix);

    gpl << cmd.str();

    std::ofstream out(prefix + ".gpl");
    out << cmd.str();
}

template<class T>
void QuantileStats<T>::plotHistAndCDF(const std::string prefix, const double bin_size, const bool absolute) const {
    plotHist(prefix + "-hist", bin_size, absolute);
    plotCDF(prefix + "-cdf");
}

template<class T>
void QuantileStats<T>::plotHistAndCDF(const std::string prefix, const double bin_size, HistConfig conf) const {
    plotHist(prefix + "-hist", bin_size, conf);
    plotCDF(prefix + "-cdf", conf);
}

template<class T>
void QuantileStats<T>::plotHistAndCDF(const std::string prefix, HistConfig conf) const {
    plotHistAndCDF(prefix, FreedmanDiaconisBinSize(), conf);
}

template<class T>
double QuantileStats<T>::FreedmanDiaconisBinSize() const {
    for (double const ignore : {1.0/2, 1.0/4, 1.0/8, 1.0/16, 1.0/32, 1.0/64}) {
        double const iqr = getQuantile(1.0 - ignore/2) - getQuantile(ignore/2);
        if (iqr > 0) {
            return 2 * iqr / cbrt(double(n));
        }
    }
    return 1;
}

template<>
double QuantileStats<size_t>::FreedmanDiaconisBinSize() const {
    for (double const ignore : {.5, .25, .125, .0625}) {
        double const iqr = getQuantile(1.0 - ignore/2) - getQuantile(.5 - ignore/2);
        if (iqr > 0) {
            return std::ceil(2 * iqr / cbrt(double(n)));
        }
    }
    return 1;
}

template<class T>
double QuantileStats<T>::getAccurateVariance() const {
    if (n < 2) {
        return 0;
    }
    double square_sum = 0;
    const double _mean = getMean();
    for (const T val : values) {
        square_sum += (val - _mean) * (val - _mean);
    }
    return square_sum / (n-1);
}

template<class T>
double QuantileStats<T>::getAccurateStddev() const {
    return std::sqrt(getAccurateVariance());
}

template<class T>
void QuantileStats<T>::reserve(const size_t size) {
    values.reserve(size);
}

template<class T>
std::vector<T> QuantileStats<T>::getData() {
    return values;
}

template<class T>
void QuantileStats<T>::clear() {
    values.clear();
    RunningStats::clear();
}

template<class T>
double QuantileStats<T>::getTrimmedMean(const T & ignore) const {
    if (ignore <= 0) {
        return getMean();
    }
    if (ignore >= 1) {
        return double(getMedian());
    }
    if (values.size() == 0) {
        return 0;
    }
    if (values.size() == 1) {
        return double(values[0]);
    }
    const size_t lower = size_t(ignore/2 * (values.size()-1));
    const size_t upper = size_t((1.0-ignore/2) * (values.size()-1));
    if (0 == upper - lower) {
        return double(getQuantile(0.5));
    }
    sort();
    double result = 0;
    for (size_t ii = lower; ii <= upper; ++ii) {
        result += double(values[ii]);
    }
    return result / (upper-lower);
}

template<class T>
template<class U>
void QuantileStats<T>::push(const std::vector<U> &_values) {
    std::lock_guard guard(push_mutex);
    push_unsafe(_values);
}

template<class T>
template<class U>
void QuantileStats<T>::push_unsafe(const std::vector<U> &_values) {
    for (const U value : _values) {
        push_unsafe(value);
    }
}

template<class T>
template<class U>
QuantileStats<U> StatsN<T>::getStats(size_t ii) const {
    if (ii >= dim()) {
        throw std::runtime_error(std::string("Requested ii (") + std::to_string(ii) + ") too large");
    }
    QuantileStats<U> result;
    result.push_unsafe(values[ii]);
    return result;
}

template QuantileStats<float> StatsN<float>::getStats(size_t) const;

template class QuantileStats<double>;
template class QuantileStats<float>;
template class QuantileStats<size_t>;

template void QuantileStats<double>::push(std::vector<double> const&);
template void QuantileStats<float>::push(std::vector<double> const&);
template void QuantileStats<size_t>::push(std::vector<double> const&);

template void QuantileStats<double>::push(std::vector<float> const&);
template void QuantileStats<float>::push(std::vector<float> const&);
template void QuantileStats<size_t>::push(std::vector<float> const&);

template void QuantileStats<double>::push_unsafe(std::vector<double> const&);
template void QuantileStats<float>::push_unsafe(std::vector<double> const&);
template void QuantileStats<size_t>::push_unsafe(std::vector<double> const&);

template void QuantileStats<double>::push_unsafe(std::vector<float> const&);
template void QuantileStats<float>::push_unsafe(std::vector<float> const&);
template void QuantileStats<size_t>::push_unsafe(std::vector<float> const&);

} // namespace runningstats
