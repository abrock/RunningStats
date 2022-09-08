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

size_t RunningStats::getCount() const {
    return n;
}

double RunningStats::getMin() const {
    return min;
}

double RunningStats::getMax() const {
    return max;
}

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

} // namespace runningstats
