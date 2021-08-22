#undef DNDEBUG
#include <cassert>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

#include "runningstats.h"

#include "gnuplot-iostream.h"

namespace runningstats {

#define READ_BIN(in, data) (in.read(reinterpret_cast<char*>(&data), sizeof data))
#define WRITE_BIN(out, data) {out.write(reinterpret_cast<const char*>(&data), sizeof data);}

RunningStats::RunningStats() {}

RunningStats::RunningStats(const RunningStats &rhs) :
    sum(rhs.sum),
    squaresum(rhs.squaresum),
    min(rhs.min),
    max(rhs.max),
    n(rhs.n),
    mean(rhs.mean),
    varSum(rhs.varSum),
    calcLog(rhs.calcLog),
    log_sum(rhs.log_sum),
    log_square_sum(rhs.log_square_sum),
    log_min(rhs.log_min),
    log_max(rhs.log_max),
    log_n(rhs.log_n)
{

}

RunningStats &RunningStats::operator=(const RunningStats &other) {
    sum = other.sum;
    squaresum = other.squaresum;
    min = other.min;
    max = other.max;
    n = other.n;
    mean = other.mean;
    varSum = other.varSum;
    calcLog = other.calcLog;
    log_sum = other.log_sum;
    log_square_sum = other.log_square_sum;
    log_min = other.log_min;
    log_max = other.log_max;
    log_n = other.log_n;
    return *this;
}

void RunningStats::clear() {
    sum = 0;
    squaresum = 0;
    n = 0;
    varSum = 0;
}

void RunningStats::push(const double value) {
    if (!std::isfinite(value)) {
        return;
    }
    const std::lock_guard<std::mutex> lock(push_mutex);
    push_unsafe(value);
}

void RunningStats::push(const std::vector<double> &data) {
    std::lock_guard guard(push_mutex);
    for (auto const val : data) {
        push_unsafe(val);
    }
}

void RunningStats::push_unsafe(const double value) {
    if (!std::isfinite(value)) {
        return;
    }
    n++;
    if (n <= 1) {
        min = value;
        max = value;
        mean = value;
        varSum = 0;
        n = 1;
    }
    else {
        min = std::min(min, double(value));
        max = std::max(max, double(value));

        sum_type const newMean = mean + (value - mean) / n;

        varSum += (value - mean) * (value - newMean);
        mean = newMean;
    }
    sum += value;
    squaresum += value * value;

    if (calcLog && value > 0) {
        double log_val = std::log(static_cast<double>(value)) / std::log(10);
        log_sum += log_val;
        log_square_sum += log_val * log_val;
        if (log_n == 0) {
            log_min = log_val;
            log_max = log_val;
        }
        else {
            log_min = std::min(log_min, log_val);
            log_max = std::max(log_max, log_val);
        }
        log_n++;
    }
}

double RunningStats::getMean() const {
    return mean;
}
double RunningStats::getLogMean() const {
    if (log_n < 1) {
        return 0;
    }
    return log_sum / log_n;
}
double RunningStats::getVar() const {
    if (n < 2) {
        return 0;
    }
    return varSum/(n - 1);
}
double RunningStats::getLogVar() const {
    if (log_n < 2) {
        return 0;
    }
    return 1.0/(log_n-1) * (log_square_sum - log_sum*log_sum / log_n);
}
double RunningStats::getStddev() const {
    return std::sqrt(getVar());
}
double RunningStats::getLogStddev() const {
    return std::sqrt(getLogVar());
}
void RunningStats::print(std::ostream& out) const {
    out.precision(15);
    out << std::scientific << getMean() << " +- " << getStddev() << ", " << n << " Samples, range: [" << min << ", " << max << "]";
}

std::string RunningStats::print() const {
    std::stringstream out;
    print(out);
    return out.str();
}

void RunningStats::printLog(std::ostream& out) const {
    out.precision(15);
    out << std::scientific << getLogMean() << " +- " << getLogStddev() << ", " << n << " Samples, range: [" << log_min << ", " << log_max << "]";
}

std::string RunningStats::printLog() const {
    std::stringstream out;
    printLog(out);
    return out.str();
}

std::string RunningStats::printBoth() const {
    return print() + "\nLogarithmic: " + printLog();
}

template<class StatsVec>
std::vector<double> RunningStats::getMean(const StatsVec& vec) {
    std::vector<double> result;
    result.reserve(vec.size());
    for (const auto& it : vec) {
        result.push_back(it.getMean());
    }
    return result;
}

template<class StatsVec>
std::vector<double> RunningStats::getStddev(const StatsVec& vec) {
    std::vector<double> result;
    result.reserve(vec.size());
    for (const auto& it : vec) {
        result.push_back(it.getStddev());
    }
    return result;
}

template<class StatsVec>
std::vector<double> RunningStats::getMin(const StatsVec& vec) {
    std::vector<double> result;
    result.reserve(vec.size());
    for (const auto& it : vec) {
        result.push_back(it.getMin());
    }
    return result;
}

template<class StatsVec>
std::vector<double> RunningStats::getMax(const StatsVec& vec) {
    std::vector<double> result;
    result.reserve(vec.size());
    for (const auto& it : vec) {
        result.push_back(it.getMax());
    }
    return result;
}

template<class StatsVec>
std::vector<size_t> RunningStats::getCount(const StatsVec& vec) {
    std::vector<size_t> result;
    result.reserve(vec.size());
    for (const auto& it : vec) {
        result.push_back(it.getCount());
    }
    return result;
}

template<class T, class StatsVec>
std::vector<T> RunningStats::getMedian(const StatsVec& vec) {
    std::vector<T> result;
    result.reserve(vec.size());
    for (const auto& it : vec) {
        result.push_back(it.getMedian());
    }
    return result;
}

template std::vector<double> RunningStats::getStddev(const std::vector<RunningStats>& vec);
template std::vector<double> RunningStats::getMean(const std::vector<RunningStats>& vec);
template std::vector<double> RunningStats::getMin(const std::vector<RunningStats>& vec);
template std::vector<double> RunningStats::getMax(const std::vector<RunningStats>& vec);
template std::vector<size_t> RunningStats::getCount(const std::vector<RunningStats>& vec);

template std::vector<double> RunningStats::getStddev(const std::vector<QuantileStats<double> >& vec);
template std::vector<double> RunningStats::getMean(const std::vector<QuantileStats<double> >& vec);
template std::vector<double> RunningStats::getMin(const std::vector<QuantileStats<double> >& vec);
template std::vector<double> RunningStats::getMax(const std::vector<QuantileStats<double> >& vec);
template std::vector<double> RunningStats::getMedian(const std::vector<QuantileStats<double> >& vec);
template std::vector<size_t> RunningStats::getCount(const std::vector<QuantileStats<double> >& vec);

size_t RunningStats::getCount() const {
    return n;
}

double RunningStats::getMin() const {
    return min;
}

double RunningStats::getMax() const {
    return max;
}

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

    return std::string("set xrange[") + std::to_string(_min) + ":" + std::to_string(_max) + "];\n";
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
    double const x_range = getMax() - getMin();
    double const n_x = x_range / _bin_size;
    if (std::isfinite(n_x) && conf.max_nx > 0 && conf.max_nx < n_x) {
        _bin_size = x_range / conf.max_nx;
    }
    Histogram h(_bin_size);
    h.push_vector_unsafe(values);
    setRangeByIgnoreAmount(conf);
    h.plotHist(prefix, conf);
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
    std::string data_file = prefix + ".data";
    cmd << "#!/usr/bin/gnuplot \n";
    cmd << "set term svg enhanced background rgb \"white\";\n"
        << "set output \"" << prefix + ".svg\"; \n"
        << conf.toString();

    if (!conf.dataLabel.empty()) {
        cmd << "set xlabel \"" << conf.dataLabel << "\";\n";
    }
    cmd << "set ylabel \"Estimated CDF\";\n";


    cmd << getXrange(conf);

    {
        std::lock_guard guard(push_mutex);
        cmd << "plot " << gpl.file(values, data_file) << " u 1:"
            << (conf.absolute ? "0" : "($0/" + std::to_string(values.size()-1) + ")") << " w l notitle; \n";
    }
    //cmd << "set term tikz; \n"
    //<< "set output \"" << prefix << ".tex\"; \n"
    //<< "replot;\n";

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
    std::string data_file = prefix + ".data";
    cmd << "#!/usr/bin/gnuplot \n";
    cmd << "set term svg enhanced background rgb \"white\";\n"
        << "set output \"" << prefix + ".svg\"; \n"
        << conf.toString();

    if (!conf.dataLabel.empty()) {
        cmd << "set xlabel \"" << conf.dataLabel << "\";\n";
    }
    cmd << "set ylabel \"Estimated CDF\";\n";

    cmd << getXrange(conf);

    cmd << "plot " << gpl.file(plot_values, data_file) << " u 1:2 w l notitle; \n";
    //cmd << "set term tikz; \n"
    // << "set output \"" << prefix << ".tex\"; \n"
    // << "replot;\n";

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
double QuantileStats<T>::FreedmanDiaconisBinSize() const {
    double const iqr = getQuantile(.75) - getQuantile(.25);
    if (iqr > 0) {
        return 2 * iqr / cbrt(double(n));
    }
    return 1;
}

template<class T>
double QuantileStats<T>::getAccurateVariance() const {
    if (n < 2) {
        return 0;
    }
    double square_sum = 0;
    const double mean = getMean();
    for (const T val : values) {
        square_sum += (val - mean) * (val - mean);
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


template class QuantileStats<double>;
template class QuantileStats<float>;

template void QuantileStats<double>::push(std::vector<double> const&);
template void QuantileStats<float>::push(std::vector<double> const&);

template void QuantileStats<double>::push(std::vector<float> const&);
template void QuantileStats<float>::push(std::vector<float> const&);

template void QuantileStats<double>::push_unsafe(std::vector<double> const&);
template void QuantileStats<float>::push_unsafe(std::vector<double> const&);

template void QuantileStats<double>::push_unsafe(std::vector<float> const&);
template void QuantileStats<float>::push_unsafe(std::vector<float> const&);

void RunningCovariance::push_unsafe(double x, double y)
{
    n++;
    if (n == 1) {
        meanX = maxX = minX = x;
        meanY = maxY = minY = y;
    }
    else {
        double const newMeanX = meanX + (x - meanX) / n;
        double const newMeanY = meanY + (y - meanY) / n;

        varSumX = varSumX + (x - meanX) * (x - newMeanX);
        varSumY = varSumY + (y - meanY) * (y - newMeanY);

        double const newCovarSum = covarSum + (x - meanX) * (y - meanY) * (n-1) / n;

        meanX = newMeanX;
        meanY = newMeanY;
        covarSum = newCovarSum;

        minX = std::min(minX, x);
        minY = std::min(minY, y);
        maxX = std::max(maxX, x);
        maxY = std::max(maxY, y);
    }
}

RunningCovariance::RunningCovariance() {}

RunningCovariance::RunningCovariance(const RunningCovariance &rhs) {
    *this = rhs;
}

RunningCovariance &RunningCovariance::operator=(const RunningCovariance &other) {
    n = other.n;
    minX = other.minX;
    maxX = other.maxX;
    minY = other.minY;
    maxY = other.maxY;
    varSumX = other.varSumX;
    varSumY = other.varSumY;
    covarSum = other.covarSum;
    return *this;
}

void RunningCovariance::push(double x, double y) {
    std::lock_guard<std::mutex> const guard(push_mutex);
    push_unsafe(x, y);
}

double RunningCovariance::getMeanX() {
    return meanX;
}

double RunningCovariance::getVarX() const {
    return n <= 1 ? 0 : varSumX/(n - 1);
}

double RunningCovariance::getMeanY() {
    return meanY;
}

double RunningCovariance::getVarY() const {
    return n <= 1 ? 0 : varSumY/(n - 1);
}

double RunningCovariance::getStddevX() const {
    return std::sqrt(getVarX());
}

double RunningCovariance::getStddevY() const {
    return std::sqrt(getVarY());
}

double RunningCovariance::getCoVar() {
    return covarSum / (n-1);
}

double RunningCovariance::getCorr() {
    return getCoVar() / std::sqrt(getVarX() * getVarY());
}

double RunningCovariance::getMinX() const {
    return minX;
}

double RunningCovariance::getMaxX() const {
    return maxX;
}

double RunningCovariance::getMinY() const {
    return minY;
}

double RunningCovariance::getMaxY() const {
    return maxY;
}

void RunningCovariance::printInfo() {
    std::cout.precision(15);
    std::cout << std::scientific
              << "Variables: " << std::endl
              << "meanX: " << meanX << std::endl
              << "meanY: " << meanY << std::endl
              << "varSumX: " << varSumX << std::endl
              << "varSumY: " << varSumY << std::endl
              << "covarSum: " << covarSum << std::endl
              << "getVarX(): " << getVarX() << std::endl
              << "getVarY(): " << getVarY() << std::endl
              << "getCoVar(): " << getCoVar() << std::endl
              << "getCorr(): " << getCorr() << std::endl
              << std::endl;
}

void BinaryStats::push(const bool value) {
    if (value) {
        count_true++;
    }
    count_total++;
}

double BinaryStats::get() const {
    if (count_total == 0) {
        return 0;
    }
    return static_cast<double>(count_true) / count_total;
}

double BinaryStats::getPercent() const {
    return 100.0*get();
}

double BinaryStats::getFalsePercent() const {
    return 100.0*(1.0-get());
}

size_t BinaryStats::getTrueCount() const {
    return count_true;
}

size_t BinaryStats::getFalseCount() const {
    return count_total - count_true;
}

std::string BinaryStats::print() const {
    std::stringstream out;
    out << "True: " << count_true << " (" << getPercent() << "%), "
        << "False: " << getFalseCount() << " (" << getFalsePercent() << "%), "
        << "Total: " << getTotalCount();
    return out.str();
}

size_t BinaryStats::getTotalCount() const {
    return count_total;
}

void BinaryStats::pushTrue(size_t const num) {
    count_true += num;
    count_total += num;
}

void BinaryStats::pushFalse(size_t const num) {
    count_total += num;
}

BinaryStats BinaryStats::merged(BinaryStats const& a, BinaryStats const& b) {
    BinaryStats result(a);
    result.pushTrue(b.getTrueCount());
    result.pushFalse(b.getFalseCount());

    return result;
}

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


bool Histogram::push_vector(const std::vector<double> &values) {
    std::lock_guard guard(push_mutex);
    return push_vector_unsafe(values);
}

bool Histogram::push_vector(const std::vector<float> &values) {
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

bool Histogram::push_vector_unsafe(const std::vector<double> &values) {
    bool success = true;
    for (const double  val : values) {
        success &= push_unsafe(val);
    }
    return success;
}

bool Histogram::push_vector_unsafe(const std::vector<float> &values) {
    bool success = true;
    for (const float  val : values) {
        success &= push_unsafe(val);
    }
    return success;
}

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
    std::string data_file = prefix + ".data";

    cmd << "#!/usr/bin/gnuplot \n";
    cmd << "set term svg enhanced background rgb \"white\";\n";
    cmd << "set output \"" << prefix + ".svg\"; \n";
    cmd << conf.toString();
    cmd << "set title \"" << conf.title << " n=" << getCount() << ", m=" << getMean() << ", s=" << getStddev() << "\"; \n";
    if (!conf.dataLabel.empty()) {
        cmd << "set xlabel \"" << conf.dataLabel << "\";\n";
    }
    cmd << "set ylabel \"Estimated PDF\";\n";

    cmd << getXrange(conf);

    auto const data = conf.absolute ? getAbsoluteHist() : getRelativeHist();

    cmd << "plot " << gpl.file1d(data, data_file) << " w boxes notitle; \n";
    //cmd << "set term tikz; \n";
    //cmd << "set output \"" << prefix << ".tex\"; \n";
    //cmd << "replot;\n";

    gpl << cmd.str();

    std::ofstream out(prefix + ".gpl");
    out << cmd.str();
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
    return std::string("set xrange[") + std::to_string(_min) + ":" + std::to_string(_max) + "];\n";
}

template<class T>
template<class U>
void QuantileStats<T>::push(const std::vector<U> &values) {
    std::lock_guard guard(push_mutex);
    push_unsafe(values);
}

template<class T>
template<class U>
void QuantileStats<T>::push_unsafe(const std::vector<U> &values) {
    for (const U value : values) {
        push_unsafe(value);
    }
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
    std::ofstream data_out(data_file);
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
    cmd << "#!/usr/bin/gnuplot \n";
    cmd << "set term svg enhanced background rgb 'white';\n";
    cmd << "set output '" << prefix << ".svg';\n";
    cmd << "plot '" << data_file << "' u 2:1:3 with image notitle;\n";
    cmd << "set term png;\n";
    cmd << "set output '" << prefix << ".png';\n";
    cmd << "replot;\n";
    gnuplotio::Gnuplot plt;
    plt << cmd.str();
    std::ofstream cmd_out(prefix + ".gpl");
    cmd_out << cmd.str();
}

template<class T>
Stats2D<T>::Stats2D(const Stats2D<T>& other) : values(other.values), quantiles_1(other.quantiles_1), quantiles_2(other.quantiles_2) {}

template<class T>
Stats2D<T>::Stats2D() {}

template<class T>
bool Stats2D<T>::push(const double a, const double b) {
    if (!std::isfinite(a) || !std::isfinite(b)) {
        return false;
    }
    std::lock_guard<std::mutex> const lock(push_mutex);
    push_unsafe(a, b);
    return true;
}

template<class T>
bool Stats2D<T>::push(const std::pair<double, double> val) {
    return push(val.first, val.second);
}

template<class T>
bool Stats2D<T>::push(const std::vector<std::pair<double, double> > &vec) {
    std::lock_guard<std::mutex> guard(push_mutex);
    return push_unsafe(vec);
}

template<class T>
bool Stats2D<T>::push_unsafe(const double a, const double b) {
    if (!std::isfinite(a) || !std::isfinite(b)) {
        return false;
    }
    values.push_back({a,b});
    quantiles_1.push_unsafe(a);
    quantiles_2.push_unsafe(b);
    return true;
}

template<class T>
bool Stats2D<T>::push_unsafe(const std::pair<double, double> val) {
    return push_unsafe(val.first, val.second);
}

template<class T>
bool Stats2D<T>::push_unsafe(const std::vector<std::pair<double, double> > &vec) {
    bool success = true;
    for (const std::pair<double, double>& val : vec) {
        success &= push_unsafe(val.first, val.second);
    }
    return success;
}

template<class T>
Histogram2D Stats2D<T>::getHistogram2D(const std::pair<double, double> bin_sizes) const {
    Histogram2D result(bin_sizes.first, bin_sizes.second);
    for (std::pair<T, T> const& d : values) {
        result.push_unsafe(d.first, d.second);
    }
    return result;
}

template<class T>
Histogram2Dfixed Stats2D<T>::getHistogram2Dfixed(const std::pair<double, double> bin_sizes, HistConfig conf) const {
    double _min_x = std::max(conf.min_x, quantiles_1.getMin());
    double _max_x = std::min(conf.max_x, quantiles_1.getMax());
    double _min_y = std::max(conf.min_y, quantiles_2.getMin());
    double _max_y = std::min(conf.max_y, quantiles_2.getMax());

    if (conf.ignore_amount > 0) {
        _min_x = std::max(_min_x, double(quantiles_1.getQuantile(conf.ignore_amount/2)));
        _max_x = std::min(_max_x, double(quantiles_1.getQuantile(1.0 - conf.ignore_amount/2)));
        _min_y = std::max(_min_y, double(quantiles_2.getQuantile(conf.ignore_amount/2)));
        _max_y = std::min(_max_y, double(quantiles_2.getQuantile(1.0 - conf.ignore_amount/2)));
    }
    Histogram2Dfixed result(
                bin_sizes.first, bin_sizes.second,
                _min_x, _min_y,
                _max_x, _max_y);
    size_t num_ignore = 0;
    for (std::pair<T, T> const& d : values) {
        bool const success = result.push_unsafe(d.first, d.second);
        if (!success) {
            num_ignore++;
        }
    }
    return result;
}

template<class T>
void Stats2D<T>::reserve(const size_t size) {
    values.reserve(size);
}

template<class T>
std::vector<std::pair<T, T> > Stats2D<T>::getData() const {
    return values;
}

template<class T>
void Stats2D<T>::plotHist(const std::string prefix, std::pair<double, double> bin_size, const HistConfig &conf) const {
    if (bin_size.first <= 0 || bin_size.second <= 0) {
        bin_size = FreedmanDiaconisBinSize();
    }
    double bin_x = bin_size.first;
    double const range_x = (quantiles_1.getMax() - quantiles_1.getMin());
    double const n_x = range_x / bin_x;
    if (conf.max_nx > 0 && n_x > conf.max_nx) {
        bin_x = range_x / conf.max_nx;
    }

    double bin_y = bin_size.second;
    double const range_y = (quantiles_2.getMax() - quantiles_2.getMin());
    double const n_y = range_y / bin_y;
    if (conf.max_ny > 0 && n_y > conf.max_ny) {
        bin_y = range_y / conf.max_ny;
    }

    getHistogram2Dfixed({bin_x, bin_y}, conf).plotHistPm3D(prefix, conf);
}

template<class T>
void Stats2D<T>::saveSummary(const std::string &filename) {
    std::ofstream out(filename);
    getSummary(out);
}

template<class T>
std::string Stats2D<T>::getSummary() {
    std::stringstream out;
    getSummary(out);
    return out.str();
}

template<class T>
void Stats2D<T>::getSummary(std::ostream &out) {
    out << "x:\n";
    quantiles_1.getSummary(out);
    out << "\ny:\n";
    quantiles_2.getSummary(out);
}

template<class T>
std::pair<double, double> Stats2D<T>::FreedmanDiaconisBinSize() const {
    double const ignore_amount = 50.0/100.0;
    double const freed_1 = quantiles_1.FreedmanDiaconisBinSize();
    double const freed_2 = quantiles_2.FreedmanDiaconisBinSize();
    double const range_1 = quantiles_1.getQuantile(1.0 - ignore_amount/2) - quantiles_1.getQuantile(ignore_amount/2);
    double const range_2 = quantiles_2.getQuantile(1.0 - ignore_amount/2) - quantiles_2.getQuantile(ignore_amount/2);
    double const num_bins_1 = range_1 / freed_1;
    double const num_bins_2 = range_2 / freed_2;
    double const num_bins = std::min(num_bins_1, num_bins_2);
    //double const num_bins_per_dim = std::sqrt(num_bins);
    double const alpha = std::sqrt(double(num_bins_1 * num_bins_2)/double(num_bins));
    std::pair<double, double> result {
        (freed_1 * alpha),
                (freed_2 * alpha)};

    /*
    double const bin_factor = std::sqrt(freed_2 / freed_1);
    std::pair<double, double> result {
        (range_1*bin_factor/num_bins_per_dim),
        (range_2/(num_bins_per_dim * bin_factor))};
        */
    return result;
}

template<class T>
size_t Stats2D<T>::size() const {
    return values.size();
}

template<class T>
bool Stats2D<T>::empty() const {
    return values.empty();
}

template<class T>
QuantileStats<T> Stats2D<T>::get1() const {
    return quantiles_1;
}

template<class T>
QuantileStats<T> Stats2D<T>::get2() const {
    return quantiles_2;
}

template<class T>
QuantileStats<T> Stats2D<T>::get(std::string const name) const {
    if (std::string("x") == name || std::string("1") == name) {
        return quantiles_1;
    }
    if (std::string("y") == name || std::string("2") == name) {
        return quantiles_2;
    }
    throw std::runtime_error("Please specify x, y, 1 or 2 in Stats2D<T>::get(std::string const name)");
}

template<class T>
std::tuple<double, double, double, double> Stats2D<T>::getQuantileEllipse(const double ignore) const {
    return {
        quantiles_1.getMedian(),
                quantiles_2.getMedian(),
                quantiles_1.getQuantile(1.0 - ignore/2) - quantiles_1.getQuantile(ignore/2),
                quantiles_2.getQuantile(1.0 - ignore/2) - quantiles_2.getQuantile(ignore/2)
    };
}

template<class T>
void Stats2D<T>::getQuantileEllipse(Ellipses &ellipse, const double ignore) const {
    ellipse.push(getQuantileEllipse(ignore));
}

template class Stats2D<float>;
template class Stats2D<double>;

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
    cmd << "set term svg enhanced background rgb 'white';\n";
    cmd << "set output '" << prefix << ".svg';\n";
    cmd << conf.toString();
    cmd << "set xrange[" << min_1 - (range_1_empty ? 1:width_1/2) << " : " << max_1 + (range_1_empty ? 1:width_1/2) << "];\n";
    cmd << "set yrange[" << min_2 - (range_2_empty ? 1:width_2/2) << " : " << max_2 + (range_2_empty ? 1:width_2/2) << "];\n";
    cmd << "set xtics out;\n";
    cmd << "set ytics out;\n";
    cmd << "plot '" << data_file << "' u 1:2:3 with image notitle;\n";
    cmd << "set term png;\n";
    cmd << "set output '" << prefix << ".png';\n";
    cmd << "replot;\n";
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
    cmd << "set term svg enhanced background rgb 'white';\n";
    cmd << "set output '" << prefix << ".svg';\n";
    cmd << conf.toString();
    cmd << "set xrange[" << min_1 - (range_1_empty ? 1:width_1/2) << " : " << max_1 + (range_1_empty ? 1:width_1/2) << "];\n";
    cmd << "set yrange[" << min_2 - (range_2_empty ? 1:width_2/2) << " : " << max_2 + (range_2_empty ? 1:width_2/2) << "];\n";
    cmd << "set key outside horiz;\n";
    cmd << "set xtics out;\n";
    cmd << "set ytics out;\n";
    cmd << "set view map;\n";
    cmd << "set pm3d corners2color c1;\n";
    cmd << "splot '" << data_file << "' u 1:2:3 with pm3d notitle";
    for (std::pair<HistConfig::Extract, double> const& extractor : conf.extractors) {
        std::vector<std::pair<double, double> > data_out;
        data_out.reserve(data.size());
        for (size_t col = 0; col < data.size(); ++col) {
            QuantileStats<float> const& col_data = stats_per_column[col];
            double const row_bin = width_1 * col + min_1;
            data_out.push_back({row_bin, col_data.getStat(extractor.first, extractor.second)});
        }
        std::string extract_name = HistConfig::extractName(extractor);
        cmd << ", " << plt.file(data_out, prefix + extract_name + ".data") << " u 1:2:(0.0) w l title '" << extract_name << "'";
    }
    cmd << ";\n";
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

HistConfig &HistConfig::setMinMaxX(const double min, const double max) {
    min_x = min;
    max_x = max;
    return *this;
}

HistConfig &HistConfig::setMinMaxY(const double min, const double max) {
    min_y = min;
    max_y = max;
    return *this;
}

HistConfig& HistConfig::setMaxPlotPts(int64_t val) {
    max_plot_pts = val;
    return *this;
}

HistConfig &HistConfig::setIgnoreAmount(const double val) {
    ignore_amount = std::min(1.0, std::max(0.0, val));
    return *this;
}

std::string HistConfig::toString() const {
    std::stringstream out;
    if (logCB) {
        out << "set logscale cb;\n";
    }
    if (logX) {
        out << "set logscale x;\n";
    }
    if (logY) {
        out << "set logscale y;\n";
    }
    if (!xLabel.empty()) {
        out << "set xlabel '" << xLabel << "';\n";
    }
    if (!yLabel.empty()) {
        out << "set ylabel '" << yLabel << "';\n";
    }
    if (!title.empty()) {
        out << "set title '" << title << "';\n";
    }
    if (fixedRatio) {
        out << "set size ratio -1;\n";
    }
    out << misc;
    return out.str();
}

std::string HistConfig::extractName(const HistConfig::Extract e) {
    switch (e) {
    case Extract::Mean: return "Mean";
    case Extract::Median: return "Median";
    case Extract::TrimmedMean: return "TrimmedMean";
    case Extract::Stddev: return "Stddev";
    case Extract::Variance: return "Variance";
    case Extract::Quantile: return "Quantile";
    case Extract::MeanAndStddev: return "MeanAndStddev";
    case Extract::MedianAndIQR: return "MedianAndIQR";
    }
    return "Unknown";
}

std::string HistConfig::extractName(const std::pair<HistConfig::Extract, double> e) {
    switch (e.first) {
    case Extract::Mean : return extractName(e.first);
    case Extract::Median : return extractName(e.first);
    case Extract::Stddev : return extractName(e.first);
    case Extract::Variance : return extractName(e.first);
    default: break;
    }
    std::string result = extractName(e.first) + "-" + std::to_string(e.second);
    size_t pos = result.size()-1;
    while (result[pos-1] == '0') {
        pos--;
    }
    return result.substr(0, pos);
}

HistConfig::Extract HistConfig::str2extract(std::string s) {
    transform(s.begin(), s.end(), s.begin(), ::tolower);
    if (s == "mean") return Extract::Mean;
    if (s == "median") return Extract::Median;
    if (s == "trimmedmean") return Extract::TrimmedMean;
    if (s == "stddev") return Extract::Stddev;
    if (s == "variance") return Extract::Variance;
    if (s == "quantile") return Extract::Quantile;
    if (s == "meandandstddev") return Extract::MeanAndStddev;
    if (s == "medianandiqr") return Extract::MedianAndIQR;
    throw std::runtime_error("Unknown extract string");
}

HistConfig& HistConfig::addExtractors(const std::vector<std::pair<std::string, double> > &vec) {
    for (auto const& it : vec) {
        addExtractors({{str2extract(it.first), it.second}});
    }
    return *this;
}

HistConfig &HistConfig::addExtractors(const std::vector<std::pair<HistConfig::Extract, double> > & vec) {
    extractors.insert(extractors.end(), vec.begin(), vec.end());
    return *this;
}

HistConfig &HistConfig::setExtractors(const std::vector<std::pair<HistConfig::Extract, double> > &vec) {
    extractors = vec;
    return *this;
}

HistConfig &HistConfig::setNormalizeX(bool value) {
    normalize_x = value;
    return *this;
}

HistConfig &HistConfig::setMaxBins(const int nx, const int ny) {
    max_nx = nx;
    max_ny = ny;
    return *this;
}

HistConfig &HistConfig::extractMean() {
    extract = Extract::Mean;
    return *this;
}
HistConfig &HistConfig::extractMedian() {
    extract = Extract::Median;
    return *this;
}
HistConfig &HistConfig::extractTrimmedMean(double const param) {
    extract = Extract::TrimmedMean;
    extractParam = param;
    return *this;
}
HistConfig &HistConfig::extractStddev() {
    extract = Extract::Stddev;
    return *this;
}
HistConfig &HistConfig::extractVariance() {
    extract = Extract::Variance;
    return *this;
}
HistConfig &HistConfig::extractQuantile(double const param) {
    extract = Extract::Quantile;
    extractParam = param;
    return *this;
}

HistConfig &HistConfig::extractMeanAndStddev() {
    extract = Extract::MeanAndStddev;
    return *this;
}

HistConfig &HistConfig::extractMedianAndIQR() {
    extract = Extract::MedianAndIQR;
    return *this;
}

HistConfig &HistConfig::setLogX(const bool val) {
    logX = val;
    return *this;
}

HistConfig &HistConfig::setLogY(const bool val) {
    logY = val;
    return *this;
}

HistConfig &HistConfig::setLogCB(bool const val) {
    logCB = val;
    return *this;
}

HistConfig &HistConfig::setAbsolute(const bool val) {
    absolute = val;
    return *this;
}

HistConfig &HistConfig::setRelative(bool const val) {
    absolute = !val;
    return *this;
}

HistConfig &HistConfig::setFixedRatio(const bool val) {
    fixedRatio = val;
    return *this;
}

HistConfig &HistConfig::setTitle(const std::string val) {
    title = val;
    return *this;
}

HistConfig &HistConfig::setXLabel(const std::string val) {
    xLabel = val;
    return *this;
}

HistConfig &HistConfig::setYLabel(const std::string val) {
    yLabel = val;
    return *this;
}

HistConfig &HistConfig::setDataLabel(const std::string val) {
    dataLabel = val;
    return *this;
}

HistConfig &HistConfig::addHorizontalLine(const double y, const std::string color) {
    misc += std::string("set arrow from graph 0, first ") + std::to_string(y) + " to graph 1, first " + std::to_string(y) + " nohead lc rgb \"" + color + "\" front;\n";
    return *this;
}

HistConfig &HistConfig::addVerticalLine(const double x, const std::string color) {
    misc += std::string("set arrow from ") + std::to_string(x) + ", graph 0 to " + std::to_string(x) + ", graph 1 nohead lc rgb \"" + color + "\" front;\n";
    return *this;
}


HistConfig HistConfig::clone() const {
    return HistConfig(*this);
}

template<class T>
Image1D<T>::Image1D(const double _width, Image2D<T> *_parent) :width(_width), parent(_parent) {}

template<class T>
T &Image1D<T>::operator[](double index) {
    int64_t ind = std::round(index / width);
    min_val = std::min(min_val, ind * width);
    max_val = std::max(max_val, ind * width);
    if (nullptr != parent) {
        parent->min_y = std::min(parent->min_y, min_val);
        parent->max_y = std::max(parent->max_y, max_val);
    }
    if (ind >= 0) {
        if (size_t(ind+2) > pos.size()) {
            pos.resize(ind+1);
        }
        return pos[ind];
    }
    size_t ind_neg = -ind;
    if (ind_neg + 2 > neg.size()) {
        neg.resize(ind_neg + 1);
    }
    return neg[ind_neg];
}

template<>
void Image1D<QuantileStats<float> >::data2file(std::ostream &out, const HistConfig &conf) {
    for (double xx = min_val; xx <= max_val + width/2; xx += width) {
        out << xx << "\t" ;
        switch (conf.extract) {
        case HistConfig::Extract::Mean: out << (this->operator[](xx)).getMean() << std::endl; break;
        case HistConfig::Extract::Stddev: out << (this->operator[](xx)).getStddev() << std::endl; break;
        case HistConfig::Extract::Variance: out << (this->operator[](xx)).getVar() << std::endl; break;
        case HistConfig::Extract::Median: out << (this->operator[](xx)).getMedian() << std::endl; break;
        case HistConfig::Extract::TrimmedMean: out << (this->operator[](xx)).getTrimmedMean(conf.extractParam) << std::endl; break;
        case HistConfig::Extract::Quantile: out << (this->operator[](xx)).getQuantile(conf.extractParam) << std::endl; break;
        case HistConfig::Extract::MeanAndStddev: out << (this->operator[](xx)).getMean() << "\t"
                                                                                         << (this->operator[](xx)).getStddev() << std::endl; break;
        case HistConfig::Extract::MedianAndIQR: out << (this->operator[](xx)).getMedian() << "\t"
                                                                                          << (this->operator[](xx)).getQuantile(.25) << "\t"
                                                                                                                                     << (this->operator[](xx)).getQuantile(.75) << std::endl; break;
        }
    }
}

template<>
void Image1D<RunningStats>::data2file(std::ostream &out, const HistConfig &conf) {
    for (double xx = min_val; xx <= max_val + width/2; xx += width) {
        out << xx << "\t" ;
        switch (conf.extract) {
        case HistConfig::Extract::Mean: out << (this->operator[](xx)).getMean() << std::endl; break;
        case HistConfig::Extract::Stddev: out << (this->operator[](xx)).getStddev() << std::endl; break;
        case HistConfig::Extract::Variance: out << (this->operator[](xx)).getVar() << std::endl; break;
        case HistConfig::Extract::MeanAndStddev: out << (this->operator[](xx)).getMean() << "\t"
                                                                                         << (this->operator[](xx)).getStddev() << std::endl; break;
        default: throw std::runtime_error("RunningStats does not provide the requested extractor");
        }
    }
    out << std::endl;
}

template<class T>
void Image1D<T>::data2file(std::ostream& out, const HistConfig &conf) {
    for (double xx = min_val; xx <= max_val + width/2; xx += width) {
        out << xx << "\t" << (this->operator[](xx)) << std::endl;
    }
}

template<class T>
void Image1D<T>::plot(const std::string &prefix, const HistConfig &conf) {
    std::string const data_file = prefix + ".data";
    std::ofstream data_out(data_file);
    data2file(data_out, conf);
    std::stringstream cmd;
    cmd << "#!/usr/bin/gnuplot \n";
    cmd << "set term svg enhanced background rgb 'white';\n";
    cmd << "set output '" << prefix << ".svg';\n";
    cmd << conf.toString() << "\n";
    cmd << "set xrange[" << min_val - width/2 << ":" << max_val + width/2 << "]; \n";
    cmd << "set xtics out;\n";
    cmd << "set ytics out;\n";
    cmd << "plot '" << data_file << "' u ";
    switch(conf.extract) {
    case HistConfig::Extract::MeanAndStddev: cmd << "1:2:3 with yerrorbars"; break;
    case HistConfig::Extract::MedianAndIQR: cmd << "1:2:3:4 with yerrorbars"; break;
    default: cmd << "1:2 w lp"; break;
    }
    cmd << " notitle;\n";
    cmd << "set term png;\n";
    cmd << "set output '" << prefix << ".png';\n";
    cmd << "replot;\n";
    gnuplotio::Gnuplot plt;
    plt << cmd.str();
    std::ofstream cmd_out(prefix + ".gpl");
    cmd_out << cmd.str();
}

template<class T>
Image2D<T>::Image2D(const double _width1, double const _width2) : width1(_width1), width2(_width2) {
}

template<class T>
Image1D<T> &Image2D<T>::operator[](const double index) {
    int64_t const ind = std::round(index / width1);
    min_x = std::min(min_x, ind * width1);
    max_x = std::max(max_x, ind * width1);
    if (ind >= 0) {
        while (size_t(ind+2) > pos.size()) {
            pos.push_back(Image1D<T>(width2, this));
        }
        return pos[ind];
    }
    size_t const ind_neg = -ind;
    while (ind_neg + 2 > neg.size()) {
        neg.push_back(Image1D<T>(width2, this));
    }
    return neg[ind_neg];
}

template<class T>
void Image2D<T>::plot(const std::string &prefix, const HistConfig &conf) {
    std::string const data_file = prefix + ".data";
    std::ofstream data_out(data_file);
    data2file(data_out, conf);
    std::stringstream cmd;
    cmd << "#!/usr/bin/gnuplot \n";
    cmd << "set term svg enhanced background rgb 'white';\n";
    cmd << "set output '" << prefix << ".svg';\n";
    cmd << conf.toString() << "\n";
    cmd << "set xrange[" << min_x - width1/2 << ":" << max_x + width1/2 << "]; \n";
    cmd << "set yrange[" << min_y - width2/2 << ":" << max_y + width2/2 << "]; \n";
    cmd << "set xtics out;\n";
    cmd << "set ytics out;\n";
    cmd << "plot '" << data_file << "' u 1:2:3 with image notitle;\n";
    cmd << "set term png;\n";
    cmd << "set output '" << prefix << ".png';\n";
    cmd << "replot;\n";
    gnuplotio::Gnuplot plt;
    plt << cmd.str();
    std::ofstream cmd_out(prefix + ".gpl");
    cmd_out << cmd.str();
}

template<>
void Image2D<RunningStats>::push_unsafe(const double x, const double y, const std::vector<double> &values) {
    for (auto const val : values) {
        (this->operator[](x))[y].push_unsafe(val);
    }
}

template<>
void Image2D<QuantileStats<float> >::push_unsafe(const double x, const double y, const std::vector<double> &values) {
    for (auto const val : values) {
        (this->operator[](x))[y].push_unsafe(val);
    }
}

template<>
void Image2D<std::vector<QuantileStats<float> > >::push_unsafe(const double x, const double y, const std::vector<double> &values) {
    std::vector<QuantileStats<float> > & target = (this->operator[](x))[y];
    if (target.size() < values.size()) {
        target.resize(values.size());
    }
    for (size_t ii = 0; ii < values.size(); ++ii) {
        target[ii].push_unsafe(values[ii]);
    }
}

template<class T>
void Image2D<T>::push_unsafe(const double x, const double y, const std::vector<double> &values) {
    for (auto const val : values) {
        (this->operator[](x))[y] += val;
    }
}

template<>
void Image2D<QuantileStats<float> >::data2file(std::ostream &out, const HistConfig &conf) {
    for (double xx = min_x; xx <= max_x + width1/2; xx += width1) {
        for (double yy = min_y; yy <= max_y + width2/2; yy += width2) {
            out << xx << "\t" << yy << "\t";
            switch (conf.extract) {
            case HistConfig::Extract::Mean: out << (this->operator[](xx))[yy].getMean() << std::endl; break;
            case HistConfig::Extract::Stddev: out << (this->operator[](xx))[yy].getStddev() << std::endl; break;
            case HistConfig::Extract::Variance: out << (this->operator[](xx))[yy].getVar() << std::endl; break;
            case HistConfig::Extract::Median: out << (this->operator[](xx))[yy].getMedian() << std::endl; break;
            case HistConfig::Extract::TrimmedMean: out << (this->operator[](xx))[yy].getTrimmedMean(conf.extractParam) << std::endl; break;
            case HistConfig::Extract::Quantile: out << (this->operator[](xx))[yy].getQuantile(conf.extractParam) << std::endl; break;
            default: throw std::runtime_error("Extractor not implemented for Image2D");
            }
        }
        out << std::endl;
    }
}

template<>
void Image2D<RunningStats>::data2file(std::ostream &out, const HistConfig &conf) {
    for (double xx = min_x; xx <= max_x + width1/2; xx += width1) {
        for (double yy = min_y; yy <= max_y + width2/2; yy += width2) {
            out << xx << "\t" << yy << "\t";
            switch (conf.extract) {
            case HistConfig::Extract::Mean: out << (this->operator[](xx))[yy].getMean() << std::endl; break;
            case HistConfig::Extract::Stddev: out << (this->operator[](xx))[yy].getStddev() << std::endl; break;
            case HistConfig::Extract::Variance: out << (this->operator[](xx))[yy].getVar() << std::endl; break;
            default: throw std::runtime_error("RunningStats does not provide the requested extractor");
            }
        }
        out << std::endl;
    }
}

template<>
void Image2D<std::vector<QuantileStats<float> > >::data2file(std::ostream &out, const HistConfig &conf) {
    for (double xx = min_x; xx <= max_x + width1/2; xx += width1) {
        for (double yy = min_y; yy <= max_y + width2/2; yy += width2) {
            out << xx << "\t" << yy;
            for (auto const& it : (this->operator[](xx))[yy]) {
                switch (conf.extract) {
                case HistConfig::Extract::Mean: out << "\t" << it.getMean(); break;
                case HistConfig::Extract::Stddev: out << "\t" << it.getStddev(); break;
                case HistConfig::Extract::Variance: out << "\t" << it.getVar(); break;
                case HistConfig::Extract::Median: out << "\t" << it.getMedian(); break;
                case HistConfig::Extract::Quantile: out << "\t" << it.getQuantile(conf.extractParam); break;
                case HistConfig::Extract::TrimmedMean: out << "\t" << it.getTrimmedMean(conf.extractParam); break;
                default: throw std::runtime_error("QuantileStats does not provide the requested extractor");
                }
            }
            out << std::endl;
        }
        out << std::endl;
    }
}

template<class T>
void Image2D<T>::data2file(std::ostream &out, const HistConfig &conf) {
    for (double xx = min_x; xx <= max_x + width1/2; xx += width1) {
        for (double yy = min_y; yy <= max_y + width2/2; yy += width2) {
            out << xx << "\t" << yy << "\t" << (this->operator[](xx))[yy] << std::endl;
        }
        out << std::endl;
    }
}

template class Image1D<double>;
template class Image1D<RunningStats>;
template class Image1D<QuantileStats<float> >;
template class Image2D<double>;
template class Image2D<float>;
template class Image2D<size_t>;
template class Image2D<RunningStats>;
template class Image2D<QuantileStats<float> >;
template class Image2D<std::vector<QuantileStats<float> > >;

template<class T>
StatsN<T>::StatsN(const std::vector<std::string> _names) : names(_names) {
    if (names.size() < 2) {
        throw std::runtime_error("Number of elements too small");
    }
}

template<class T>
StatsN<T>::StatsN(const StatsN<T> &other) : names(other.names), values(other.values) {}

template<class T>
StatsN<T> &StatsN<T>::operator =(const StatsN<T> &other) {
    names = other.names;
    values = other.values;
    return *this;
}

template<class T>
size_t StatsN<T>::dim() const {
    return names.size();
}

template<class T>
Stats2D<T> StatsN<T>::getStats2D(size_t ii, size_t jj) const {
    if (ii >= dim()) {
        throw std::runtime_error(std::string("Requested ii (") + std::to_string(ii) + ") too large");
    }
    if (jj >= dim()) {
        throw std::runtime_error(std::string("Requested jj (") + std::to_string(jj) + ") too large");
    }
    if (ii == jj) {
        throw std::runtime_error(std::string("ii and jj are identical (") + std::to_string(ii) + ")");
    }
    Stats2D<T> result;
    result.reserve(size());

    for (std::vector<T> const& v : values) {
        result.push_unsafe(v[ii], v[jj]);
    }

    return result;
}

template<class T>
void StatsN<T>::plotAll(const std::string prefix, const HistConfig conf) const {
    HistConfig c(conf);
    for (size_t ii = 1; ii < dim(); ++ii) {
        std::string const& x_label = names[ii];
        c.setXLabel(x_label);
        for (size_t jj = 0; jj < ii; ++jj) {
            std::string const& y_label = names[jj];
            c.setYLabel(y_label);
            Stats2D<T> stats = getStats2D(ii, jj);
            stats.plotHist(prefix + "-" + x_label + "-" + y_label, stats.FreedmanDiaconisBinSize(), c);
        }
    }
}

template<class T>
size_t StatsN<T>::size() const {
    return values.size();
}

template class StatsN<float>;

template<class T>
template<class U>
bool StatsN<T>::push_unsafe(const std::vector<U> &val) {
    if (dim() != val.size()) {
        throw std::runtime_error(std::string("Size of given vector (") + std::to_string(val.size()) + ") doesn't match StatN object size (" + std::to_string(dim()) + ")");
    }
    for (const U v: val) {
        if (!std::isfinite(v)) {
            return false;
        }
    }
    values.push_back(std::vector<T> (val.begin(), val.end()));
    return true;
}

template bool StatsN<float>::push_unsafe(const std::vector<float> &val);
template bool StatsN<float>::push_unsafe(const std::vector<double> &val);

template<class T>
template<class U>
bool StatsN<T>::push(const std::vector<U> &val) {
    std::lock_guard<std::mutex> guard(push_mutex);
    return push_unsafe(val);
}

template bool StatsN<float>::push(const std::vector<float> &val);
template bool StatsN<float>::push(const std::vector<double> &val);

template<class T>
template<class U>
bool Stats2D<T>::push_unsafe(const Stats2D<U> &other) {
    bool success = true;
    for (auto const& val : other.getData()) {
        success &= push_unsafe(val.first, val.second);
    }
    return success;
}

template<class T>
template<class U>
bool Stats2D<T>::push(const Stats2D<U> &other) {
    std::lock_guard<std::mutex> guard(push_mutex);
    return push_unsafe(other);
}

template bool Stats2D<float>::push_unsafe(const Stats2D<float> & other);
template bool Stats2D<float>::push_unsafe(const Stats2D<double> & other);
template bool Stats2D<double>::push_unsafe(const Stats2D<float> & other);
template bool Stats2D<double>::push_unsafe(const Stats2D<double> & other);

template bool Stats2D<float>::push(const Stats2D<float> & other);
template bool Stats2D<float>::push(const Stats2D<double> & other);
template bool Stats2D<double>::push(const Stats2D<float> & other);
template bool Stats2D<double>::push(const Stats2D<double> & other);

void Ellipses::push_unsafe(const std::tuple<double, double, double, double> &val) {
    data.push_back(val);

    limits_x.push_unsafe(std::get<0>(val) + std::get<2>(val)/2);
    limits_x.push_unsafe(std::get<0>(val) - std::get<2>(val)/2);

    limits_y.push_unsafe(std::get<1>(val) + std::get<3>(val)/2);
    limits_y.push_unsafe(std::get<1>(val) - std::get<3>(val)/2);
}

void Ellipses::push(const std::tuple<double, double, double, double> &val) {
    std::lock_guard<std::mutex> guard(push_mutex);
    push_unsafe(val);
}

void Ellipses::plot(std::string const& prefix, const HistConfig &conf) {
    gnuplotio::Gnuplot gpl;
    std::stringstream cmd;
    std::string data_file = prefix + ".data";
    double const size_x = limits_x.getMax() - limits_x.getMin();
    double const size_y = limits_y.getMax() - limits_y.getMin();
    double const margin = 0.04;
    cmd << "#!/usr/bin/gnuplot \n";
    cmd << "set term svg enhanced background rgb \"white\";\n"
        << "set output \"" << prefix + ".svg\"; \n"
        << conf.toString() << ";\n"
        << "set xrange [" << limits_x.getMin() - margin * size_x << ":" << limits_x.getMax() + margin * size_x << "];\n"
        << "set yrange [" << limits_y.getMin() - margin * size_y << ":" << limits_y.getMax() + margin * size_y << "];\n"
           ;

    cmd << "plot " << gpl.file(data, data_file) << " with ellipses notitle; \n";
    //cmd << "set term tikz; \n"
    //<< "set output \"" << prefix << ".tex\"; \n"
    //<< "replot;\n";

    gpl << cmd.str();

    std::ofstream out(prefix + ".gpl");
    out << cmd.str();
}

WeightedRunningStats::WeightedRunningStats() {}

WeightedRunningStats::WeightedRunningStats(const WeightedRunningStats &rhs)
    :
      sum(rhs.sum),
      squaresum(rhs.squaresum),
      weight_sum(rhs.weight_sum),
      min(rhs.min),
      max(rhs.max),
      n(rhs.n),
      mean(rhs.mean),
      varSum(rhs.varSum)
{}

void WeightedRunningStats::push_unsafe(const double val, const double weight) {
    sum += val * weight;
    squaresum += val * (val * weight);
    weight_sum += weight;
    n++;
}

double WeightedRunningStats::getMean() const {
    if (n <= 0 || weight_sum <= 0) {
        return 0;
    }
    return sum / weight_sum;
}

double WeightedRunningStats::getVar() const {
    if (n <= 0 || weight_sum <= 0) {
        return 0;
    }
    double const mean = getMean();
    return (squaresum - 2 * mean * sum + mean * mean * weight_sum) / weight_sum;
}

double WeightedRunningStats::getStddev() const {
    return std::sqrt(getVar());
}

template<class T>
void ThresholdErrorMean<T>::push_unsafe(const T val, const T error) {
    data.push_back({val, error});
}

template<class T>
void ThresholdErrorMean<T>::plot(const std::string &prefix, const HistConfig &conf) const {
    if (data.size() < 2) {
        return;
    }
    std::lock_guard guard(access_lock);
    std::vector<T> only_errors;
    only_errors.reserve(data.size());
    for (size_t ii = 0; ii < data.size(); ++ii) {
        only_errors.push_back(data[ii].second);
    }
    std::sort(data.begin(), data.end());
    std::sort(only_errors.begin(), only_errors.end());

    std::vector<std::pair<HistConfig::Extract, double> > extractors = conf.extractors;
    if (extractors.empty()) {
        extractors.push_back({conf.extract, conf.extractParam});
    }

    std::string const error_file1 = prefix + "-errors-over-threshold";
    std::string const percentage_file1 = prefix + "-percentages-over-threshold";

    gnuplotio::Gnuplot plt1;
    std::stringstream cmd1;
    cmd1 << "#!/usr/bin/gnuplot \n";
    cmd1 << "set term svg enhanced background rgb 'white';\n";
    cmd1 << "set output '" << prefix << "-threshold.svg';\n";
    cmd1 << "set key out horiz;\n";
    cmd1 << conf.toString();
    cmd1 << "set y2tics;\n set ytics nomirror;\n";
    cmd1 << "plot ";

    std::string const error_file2 = prefix + "-errors-over-percentage";
    std::string const percentage_file2 = prefix + "-thresholds-over-percentage";

    gnuplotio::Gnuplot plt2;
    std::stringstream cmd2;
    cmd2 << "#!/usr/bin/gnuplot \n";
    cmd2 << "set term svg enhanced background rgb 'white';\n";
    cmd2 << "set output '" << prefix << "-percentage.svg';\n";
    cmd2 << "set key out horiz;\n";
    cmd2 << conf.toString();
    cmd2 << "set xlabel 'Percentage of data kept';\n";
    cmd2 << "set ylabel '';\n";
    cmd2 << "set y2tics;\n set ytics nomirror;\n";
    cmd2 << "plot ";


    bool first_extract = true;

    std::vector<std::pair<T, T> >
            error_over_threshold, percentage_over_threshold,
            error_over_percentage, error_over_percentage_oracle, threshold_over_percentage;
    for (std::pair<HistConfig::Extract, double> const& extractor : extractors) {
        error_over_threshold.clear();
        percentage_over_threshold.clear();
        error_over_percentage.clear();
        threshold_over_percentage.clear();
        error_over_percentage_oracle.clear();
        QuantileStats<T> error_stats, oracle_stats;
        for (size_t ii = 0; ii < data.size(); ++ii) {
            auto const& it = data[ii];
            oracle_stats.push_unsafe(only_errors[ii]);
            error_stats.push_unsafe(it.second);
            if (conf.max_plot_pts <= 0
                    || conf.max_plot_pts >= int64_t(data.size())
                    || double(ii) * (double(conf.max_plot_pts) / data.size()) >= error_over_threshold.size()) {
                double const error = error_stats.getStat(extractor.first, extractor.second);
                double const oracle_error = oracle_stats.getStat(extractor.first, extractor.second);
                double const threshold = it.first;
                double const percentage = 100.0*double(ii) / double(data.size());
                error_over_threshold.push_back({threshold, error});
                percentage_over_threshold.push_back({threshold, percentage});
                error_over_percentage.push_back({percentage, error});
                error_over_percentage_oracle.push_back({percentage, oracle_error});
                threshold_over_percentage.push_back({percentage, threshold});
            }
        }
        double const error = error_stats.getStat(extractor.first, extractor.second);
        double const threshold = data.rbegin()->first;
        double const percentage = 100.0;
        error_over_threshold.push_back({threshold, error});
        percentage_over_threshold.push_back({threshold, percentage});
        error_over_percentage.push_back({percentage, error});
        threshold_over_percentage.push_back({percentage, threshold});
        std::string const extractName = conf.extractName(extractor);
        cmd1 << plt1.file(error_over_threshold, error_file1 + extractName + ".data")
             << " u 1:2 with l title 'error " << extractName << "', ";
        cmd2 << plt2.file(error_over_percentage, error_file2 + extractName + ".data")
             << " u 1:2 with l title 'error " << extractName << "',";
        cmd2 << plt2.file(error_over_percentage_oracle, error_file2 + extractName + "-oracle.data")
             << " u 1:2 with l title 'oracle " << extractName << "',";

        if (first_extract) {
            std::cout << "Error stats: " << error_stats.print() << std::endl;
            error_stats.plotHistAndCDF(prefix + "-errors", -1, HistConfig()
                                       .setMaxPlotPts(1000)
                                       .setMaxBins(1000,1000)
                                       .setDataLabel("Errors"));
        }
        first_extract = false;
    }

    cmd1 << plt1.file(percentage_over_threshold, percentage_file1) << " u 1:2 axes x1y2 w l title 'percentage';\n";
    cmd1 << "set term png;\n";
    cmd1 << "set output '" << prefix << "-error.png';\n";
    cmd1 << "replot;\n";
    plt1 << cmd1.str();
    std::ofstream cmd_out1(prefix + "-over-threshold.gpl");
    cmd_out1 << cmd1.str();

    cmd2 << plt2.file(threshold_over_percentage, percentage_file2) << " u 1:2 axes x1y2 w l title 'threshold';\n";
    cmd2 << "set term png;\n";
    cmd2 << "set output '" << prefix << "-percentage.png';\n";
    cmd2 << "replot;\n";
    plt2 << cmd2.str();
    std::ofstream cmd_out2(prefix + "-over-percentage.gpl");
    cmd_out2 << cmd2.str();

}

template<class T>
void ThresholdErrorMean<T>::save(std::ostream &out) const {
    std::lock_guard guard(access_lock);
    for (auto const& it : data) {
        WRITE_BIN(out, it.first);
        WRITE_BIN(out, it.second);
    }
}

template<class T>
void ThresholdErrorMean<T>::save(const std::string &filename) const {
    std::ofstream out(filename);
    save(out);
}

template<class T>
void ThresholdErrorMean<T>::load(std::istream &in) {
    std::lock_guard guard(access_lock);
    while (in) {
        T x = 0, y = 0;
        if (!READ_BIN(in, x)) {
            return;
        }
        if (!READ_BIN(in, y)) {
            return;
        }
        push_unsafe(x, y);
    }
}

template<class T>
void ThresholdErrorMean<T>::load(const std::string &filename) {
    std::ifstream in(filename);
    load(in);
}

template<class T>
size_t ThresholdErrorMean<T>::size() const {
    return data.size();
}

template<class T>
std::vector<std::pair<T, T> > ThresholdErrorMean<T>::getData() const {
    return data;
}

template class ThresholdErrorMean<float>;

} // namespace runningstats
