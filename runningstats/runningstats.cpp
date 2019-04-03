#include <algorithm>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

#include "runningstats.h"

namespace runningstats {

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
#pragma omp critical
    {
        push_unsafe(value);
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

        double const newMean = mean + (value - mean) / n;

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
    out << getMean() << " +- " << getStddev() << ", " << n << " Samples, range: [" << min << ", " << max << "]";
}

std::string RunningStats::print() const {
    std::stringstream out;
    print(out);
    return out.str();
}

void RunningStats::printLog(std::ostream& out) const {
    out << getLogMean() << " +- " << getLogStddev() << ", " << n << " Samples, range: [" << log_min << ", " << log_max << "]";
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
    sorted = false;
    values.push_back(value);
    RunningStats::push(value);
}


template<class T>
T QuantileStats<T>::getQuantile(const double quantile) const {
    if (quantile <= 0) {
        return min;
    }
    if (quantile >= 1) {
        return max;
    }
    if (values.size() == 0) {
        return 0;
    }
    if (values.size() == 1) {
        return values[0];
    }
    sort();
    return values[static_cast<size_t>(quantile * (values.size()-1))];
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
        std::sort(values.begin(), values.end());
        sorted = true;
    }
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
    sort();
    return values;
}

template<class T>
double QuantileStats<T>::getTrimmedMean(const T & ignore) {
    if (ignore <= 0) {
        return getMean();
    }
    if (ignore >= 1) {
        return double(getQuantile(0.5));
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

void RunningCovariance::push(double x, double y)
{
    n++;
    if (n == 1) {
        meanX = x;
        meanY = y;
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
    }
}

double RunningCovariance::getMeanX() {
    return meanX;
}

double RunningCovariance::getVarX() {
    return n <= 1 ? 0 : varSumX/(n - 1);
}

double RunningCovariance::getMeanY() {
    return meanY;
}

double RunningCovariance::getVarY() {
    return n <= 1 ? 0 : varSumY/(n - 1);
}

double RunningCovariance::getCoVar() {
    return covarSum / (n-1);
}

double RunningCovariance::getCorr() {
    return getCoVar() / std::sqrt(getVarX() * getVarY());
}

void RunningCovariance::printInfo()
{
    std::cout << "Variables: " << std::endl
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
    return 100*get();
}

size_t BinaryStats::getTrueCount() const {
    return count_true;
}

size_t BinaryStats::getFalseCount() const {
    return count_total - count_true;
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

bool Histogram::push(double value) {
    if (!std::isfinite(value)) {
        return false;
    }
#pragma omp critical
    {
        push_unsafe(value);
    }
    return true;
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

std::vector<std::pair<double, double> > Histogram::getAbsoluteHist() const {
    std::vector<std::pair<double, double> > result;
    for (auto it = data.begin(); it != data.end(); ++it) {
        result.push_back(std::pair<double, double>(it->first * bin_size, it->second));
    }
    return result;
}

std::vector<std::pair<double, double> > Histogram::getRelativeHist() const {
    std::vector<std::pair<double, double> > result;
    for (auto it = data.begin(); it != data.end(); ++it) {
        result.push_back(std::pair<double, double>(it->first * bin_size, double(it->second)/n));
    }
    return result;
}

} // namespace runningstats
