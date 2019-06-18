#ifndef RUNNINGSTATS_H
#define RUNNINGSTATS_H

#include <algorithm>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <map>

namespace runningstats {

/**
 * @brief The BinaryStats class collects statistics about binary values.
 */
class BinaryStats {
private:
    size_t count_total = 0;
    size_t count_true = 0;
public:
    /**
     * @brief push Adds a binary value to the statistics.
     * @param value
     */
    void push(bool const value);
    /**
     * @brief get returns the number of "true" pushes divided by the total number of pushes. Range is [0,1]
     * @return
     */
    double get() const;
    /**
     * @brief getPercent returns the percentage of "true" pushes. Range is [0,100]
     * @return
     */
    double getPercent() const;
    size_t getTrueCount() const;
    size_t getTotalCount() const;
    static BinaryStats merged(const BinaryStats &a, const BinaryStats &b);
    void pushTrue(const size_t num);
    void pushFalse(const size_t num);
    size_t getFalseCount() const;
};

class RunningCovariance
{
public:
    size_t n = 0;
    double
    minX = 0, maxX = 0,
    minY = 0, maxY = 0,
    meanX = 0, meanY = 0,
    varSumX = 0, varSumY = 0,
    covarSum = 0;

    void push(double x, double y);

    size_t getN()
    {
        return n;
    }

    double getMeanX();
    double getVarX();

    double getMeanY();
    double getVarY();

    double getCoVar();
    double getCorr();

    double getMinX() const;
    double getMaxX() const;
    double getMinY() const;
    double getMaxY() const;

    void printInfo();
};

class RunningStats {

public:

    template<class StatsVec>
    static std::vector<double> getMean(const StatsVec& vec);

    template<class StatsVec>
    static std::vector<double> getStddev(const StatsVec& vec);

    template<class StatsVec>
    static std::vector<double> getMin(const StatsVec& vec);

    template<class StatsVec>
    static std::vector<double> getMax(const StatsVec& vec);

    template<class StatsVec>
    static std::vector<size_t> getCount(const StatsVec& vec);

    template<class T, class StatsVec>
    static std::vector<T> getMedian(const StatsVec& vec);

    void clear();

    /**
     * @brief push add a value to the statistics. This basically calls push_unsafe inside a #pragma omp critical block.
     * @param value
     */
    void push(const double value);

    /**
     * @brief push_unsafe adds a value to the statistics. This function must not be called from multiple threads at the same time.
     * @param value
     */
    void push_unsafe(const double value);

    double getMean() const;

    double getLogMean() const;

    double getVar() const;

    double getLogVar() const;

    double getStddev() const;

    double getLogStddev() const;

    double getMin() const;

    double getMax() const;

    size_t getCount() const;

    void print(std::ostream& out) const;

    std::string print() const;

    void printLog(std::ostream& out) const;

    std::string printLog() const;

    std::string printBoth() const;

    double sum = 0;
    double squaresum = 0;
    double min = 0;
    double max = 0;
    size_t n = 0;

    double mean = 0;
    double varSum = 0;

    bool calcLog = true;

    double log_sum = 0;
    double log_square_sum = 0;
    double log_min = 0;
    double log_max = 0;
    size_t log_n = 0;
};


class Histogram : public RunningStats {
    double const bin_size;
    std::map<int64_t, size_t> data;

public:
    Histogram(double const _bin_size);

    bool push(double const value);

    bool push_unsafe(double const value);

    std::vector<std::pair<double, double> > getAbsoluteHist() const;
    std::vector<std::pair<double, double> > getRelativeHist() const;
};

template<class T>
class QuantileStats : public RunningStats {
public:
    void push(const double value);

    T getQuantile(const double quantile) const;

    T getInverseQuantile(const double value) const;

    T getMedian() const;

    double getAccurateVariance() const;

    double getAccurateStddev() const;

    void reserve(const size_t size);

    std::vector<T> getData();

    double getTrimmedMean(const T & ignore);

private:
    void sort() const;

    mutable std::vector<T> values;
    mutable bool sorted = true;
};

} // namespace runningstats

#endif // RUNNINGSTATS_H
