#ifndef RUNNINGSTATS_H
#define RUNNINGSTATS_H

#include <algorithm>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <map>
#include <set>

#include <mutex>

#include <boost/multiprecision/float128.hpp>

namespace runningstats {

class HistConfig {
public:
    /**
     * @brief logX logarithmic scale of the x axis.
     */
    bool logX = false;

    /**
     * @brief logY logarithmic scale of the y axis.
     */
    bool logY = false;

    /**
     * @brief logCB logarithmic scale of the heatmap (for 2D histograms only)
     */
    bool logCB = false;

    /**
     * @brief absolute Use absolute bin counts instead of estimated probability densities
     */
    bool absolute = false;

    /**
     * @brief Set fixedRatio to true to make Gnuplot plot both axes at the same scale.
     */
    bool fixedRatio = false;

    double min_x = -std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::max();

    double min_y = -std::numeric_limits<double>::max();
    double max_y = std::numeric_limits<double>::max();

    double ignore_amount = 0;

    /**
     * @brief max_nx Maximum number of bins in x-direction
     */
    int max_nx = -1;
    /**
     * @brief max_ny Maximum number of bins in y-direction
     */
    int max_ny = -1;

    /**
     * @brief max_plot_pts is the maximum number of values to be used in a plot.
     */
    int64_t max_plot_pts = -1;

    HistConfig &setMaxPlotPts(int64_t val);

    HistConfig & setIgnoreAmount(double const val);

    std::string title;

    std::string xLabel;

    std::string yLabel;

    std::string dataLabel;

    std::string toString() const;

    std::string misc;

    enum class Extract {
        Mean,
        Median,
        TrimmedMean,
        Stddev,
        Variance,
        Quantile,
        MeanAndStddev,
        MedianAndIQR
    };

    /**
     * @brief extractName returns a string representing an Extract option
     * @param e
     * @return
     */
    static std::string extractName(Extract const e);

    /**
     * @brief extractName returns a string representing an Extract option + parameter
     * @param e
     * @return
     */
    static std::string extractName(std::pair<Extract, double> const e);

    /**
     * @brief str2extract finds the Extract value for a given string or throws an expcetion.
     * @param s
     * @return
     */
    static Extract str2extract(std::string s);

    /**
     * @brief extractors Set of extractors to be plotted in ThresholdErrorMean.
     */
    std::vector<std::pair<Extract, double> > extractors;

    HistConfig& addExtractors(std::vector<std::pair<Extract, double> > const& vec);
    HistConfig& setExtractors(std::vector<std::pair<Extract, double> > const& vec);

    HistConfig& addExtractors(std::vector<std::pair<std::string, double> > const& vec);

    Extract extract = Extract::Mean;
    double extractParam = .5;

    /**
     * @brief normalize_x if this is set to true then for each bin in y-direction the values in x-direction will be normalized.
     */
    bool normalize_x = false;
    HistConfig& setNormalizeX(bool value = true);

    HistConfig& setMaxBins(int const nx, int const ny);

    HistConfig& extractMean();
    HistConfig& extractMedian();
    HistConfig& extractTrimmedMean(double const param = 0.5);
    HistConfig& extractStddev();
    HistConfig& extractVariance();
    HistConfig& extractQuantile(double const param = 0.5);
    HistConfig& extractMeanAndStddev();
    HistConfig& extractMedianAndIQR();

    HistConfig& setMinMaxX(double const min, double const max);
    HistConfig& setMinMaxY(double const min, double const max);

    HistConfig& setLogX(bool const val = true);
    HistConfig& setLogY(bool const val = true);
    HistConfig& setLogCB(bool const val = true);
    HistConfig& setAbsolute(bool const val = true);
    HistConfig& setRelative(bool const val = true);
    HistConfig& setFixedRatio(bool const val = true);

    HistConfig& setTitle(std::string const val);
    HistConfig& setXLabel(std::string const val);
    HistConfig& setYLabel(std::string const val);
    HistConfig& setDataLabel(std::string const val);

    HistConfig &addHorizontalLine(double const y, std::string const color = "#ffffff");
    HistConfig &addVerticalLine(double const x, std::string const color = "#ffffff");

    HistConfig clone() const;
};

template<class T>
/**
 * @brief The ThresholdErrorMean class collects pairs of values consisting of a confidence measure (larger is better)
 * and an error measure (smaller is better) and creates plots to answer the question:
 * If I discard pairs with a confidence lower than x, what is the mean error and what is the percentage of values dropped?
 */
class ThresholdErrorMean {
public:
    void push_unsafe(const T val, const T error);

    void plot(std::string const& prefix, HistConfig const& conf = HistConfig()) const;

    void save(std::ostream& out) const;
    void save(std::string const& filename) const;

    void load(std::istream &in);
    void load(std::string const &filename);

    size_t size() const;

    std::vector<std::pair<T, T> > getData() const;
private:
    mutable std::mutex access_lock;
    mutable std::vector<std::pair<T, T> > data;
};


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
    double getFalsePercent() const;
    size_t getTrueCount() const;
    size_t getTotalCount() const;
    static BinaryStats merged(const BinaryStats &a, const BinaryStats &b);
    void pushTrue(const size_t num);
    void pushFalse(const size_t num);
    size_t getFalseCount() const;

    std::string print() const;
};

class RunningCovariance {
    std::mutex push_mutex;
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
    double getVarX() const;

    double getMeanY();
    double getVarY() const;

    double getStddevX() const;
    double getStddevY() const;

    double getCoVar();
    double getCorr();

    double getMinX() const;
    double getMaxX() const;
    double getMinY() const;
    double getMaxY() const;

    void printInfo();
    void push_unsafe(double x, double y);
};

class WeightedRunningStats {
private:
    std::mutex push_mutex;

    double sum = 0;
    double squaresum = 0;
    double weight_sum = 0;
    double min = 0;
    double max = 0;
    size_t n = 0;

    double mean = 0;
    double varSum = 0;

public:

    WeightedRunningStats();
    WeightedRunningStats(WeightedRunningStats const& rhs);
    WeightedRunningStats & operator= ( const WeightedRunningStats & other);

    void push_unsafe(double const val, double const weight);

    double getMean() const;

    double getVar() const;

    double getStddev() const;
};

class RunningStats {

protected:
    mutable std::mutex push_mutex;

public:

    RunningStats();
    RunningStats(RunningStats const& rhs);
    RunningStats & operator= ( const RunningStats & other);

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

    void push(const std::vector<double>& data);

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

    typedef double sum_type;

    sum_type sum = 0;
    sum_type squaresum = 0;
    double min = 0;
    double max = 0;
    size_t n = 0;

    sum_type mean = 0;
    sum_type varSum = 0;

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

    Histogram(Histogram const& rhs);

    bool push(double const value);

    bool push_vector(std::vector<double> const& values);
    bool push_vector(std::vector<float> const& values);

    bool push_unsafe(double const value);

    bool push_vector_unsafe(std::vector<double> const& values);
    bool push_vector_unsafe(std::vector<float> const& values);

    std::vector<std::pair<double, double> > getAbsoluteHist() const;
    std::vector<std::pair<double, double> > getRelativeHist() const;

    void plotHist(std::string const prefix, bool const absolute = true) const;
    void plotHist(std::string const prefix, HistConfig const conf) const;

    size_t getBinCount(double const value) const;

    double getLikelihood() const;
    double getPosteriorProbability() const;

    std::string getXrange(HistConfig const& conf) const;
};

class Histogram2D {
    double const width_1;
    double const width_2;

    std::set<int64_t> bins_1, bins_2;

    /**
     * @brief data Bin counts can be addressed as data[row][col];
     */
    std::map<int64_t, std::map<int64_t, size_t> > data;

    std::mutex push_mutex;

    size_t total_count = 0;

public:
    Histogram2D(double const _bin_1, double const _bin_2);
    Histogram2D(Histogram2D const& rhs);

    bool push(double const val1, double const val2);

    bool push_unsafe(double const val1, double const val2);

    void plotHist(std::string const prefix, bool const absolute = true) const;
};

template<class T>
class QuantileStats;

class Histogram2Dfixed {
    double const width_1;
    double const width_2;
    double const min_1, min_2;
    double const max_1, max_2;

    /**
     * @brief data Bin counts can be addressed as data[row][col];
     */
    std::vector<std::vector<size_t> > data;
    std::vector<QuantileStats<float> > stats_per_column;

    std::mutex push_mutex;

    size_t total_count = 0;


public:
    Histogram2Dfixed(double const _bin_1, double const _bin_2,
                     double const _min_1, double const _min_2,
                     double const _max_1, double const _max_2);
    Histogram2Dfixed(Histogram2Dfixed const& rhs);

    ~Histogram2Dfixed();

    static double sanitize_bin_width(double const width, double const min, double const max);

    bool push(double const val1, double const val2);

    bool push_unsafe(double const val1, double const val2);

    void plotHist(std::string const prefix, const bool absolute = true) const;
    void plotHist(std::string const prefix, const HistConfig &conf = HistConfig()) const;

    void plotHistPm3D(std::string const prefix, const HistConfig &conf = HistConfig()) const;
};


template<class T>
class QuantileStats : public RunningStats {
public:
    void push(const double value);

    void push_unsafe(const double value);

    template<class U>
    void push(std::vector<U> const& values);

    template<class U>
    void push_unsafe(std::vector<U> const& values);

    Histogram getHistogram(double const bin_size);

    T getQuantile(const double quantile) const;

    T getInverseQuantile(const double value) const;

    T getMedian() const;

    double getAccurateVariance() const;

    double getAccurateStddev() const;

    void reserve(const size_t size);

    std::vector<T> getData();

    /**
     * @brief getTrimmedMean calculates the mean of the central part of the data, ignoring top and bottom outliers.
     * @param ignore Amount of data to ignore. 0.5 means ignore 50%, so the top 25% and bottom 25% is ignored.
     * @return
     */
    double getTrimmedMean(const T & ignore) const;

    void sort() const;

    void plotHist(std::string const prefix, double const bin_size, const bool absolute = true) const;
    void plotHist(std::string const prefix, double const bin_size, HistConfig conf) const;

    void plotCDF(std::string const prefix, HistConfig conf = HistConfig()) const;

    /**
     * @brief plotReducedCDF plots a CDF (cumulative distribution function) which should look exactly
     * like the one plotted by #plotCDF but doesn't use more than 5k points. This should be sufficient
     * to print the resulting plot on a DIN A2 poster at 300 dpi without noticable difference to the plot by
     * #plotCDF
     * @param prefix
     */
    void plotReducedCDF(std::string const prefix, HistConfig conf = HistConfig()) const;

    void plotHistAndCDF(std::string const prefix, double const bin_size, const bool absolute = true) const;

    void plotHistAndCDF(std::string const prefix, double const bin_size, HistConfig conf) const;

    double FreedmanDiaconisBinSize() const;

    static T getQuantile(const double quantile, std::vector<T> &values, bool &sorted);

    static double getMin(const std::vector<T> &values);
    static double getMax(const std::vector<T> &values);

    std::string getXrange(HistConfig const& conf) const;

    void setRangeByIgnoreAmount(HistConfig & conf) const;

    using RunningStats::getMin;
    using RunningStats::getMax;

    void saveSummary(std::string const & filename);

    std::string getSummary();

    void getSummary(std::ostream& out);

    double getStat(HistConfig::Extract const type, double const param) const;
    double getStat(const std::pair<HistConfig::Extract, double> pair) const;

    void save(std::ostream &out) const;
    void save(std::string const& filename) const;

    void load(std::istream &in);
    void load(std::string const & filename);

private:

    mutable std::vector<T> values;
    mutable bool sorted = true;
};

/**
 * @brief The Ellipses class stores ellipses defined by center point and x/y diameters. It plots them.
 */
class Ellipses {
private:
    std::vector<std::tuple<double, double, double, double> > data;

    RunningStats limits_x;
    RunningStats limits_y;

    std::mutex push_mutex;
public:
    void push_unsafe(std::tuple<double, double, double, double> const& val);
    void push(std::tuple<double, double, double, double> const& val);

    void plot(const std::string &prefix, HistConfig const& conf);
};

template<class T>
class Stats2D {
private:
    std::mutex push_mutex;
public:
    Stats2D(const Stats2D<T> &other);
    Stats2D();

    bool push(const double a, const double b);
    bool push(const std::pair<double, double> val);
    bool push(const std::vector<std::pair<double, double> >& vec);
    template<class U>
    bool push(const Stats2D<U> & other);

    bool push_unsafe(const double a, const double b);
    bool push_unsafe(const std::pair<double, double> val);
    bool push_unsafe(const std::vector<std::pair<double, double> >& vec);
    template<class U>
    bool push_unsafe(const Stats2D<U> & other);

    Histogram2D getHistogram2D(std::pair<double, double> const bin_sizes) const;

    Histogram2Dfixed getHistogram2Dfixed(std::pair<double, double> const bin_sizes, HistConfig conf = HistConfig()) const;

    std::pair<double, double> getMedian() const;

    void reserve(const size_t size);

    std::vector<std::pair<T, T> > getData() const;

    double getTrimmedMean(const T & ignore);

    void sort() const;

    void plotHist(std::string const prefix, std::pair<double, double> bin_size = {-1,-1}, const HistConfig &conf = HistConfig()) const;

    void saveSummary(std::string const & filename);

    std::string getSummary();

    void getSummary(std::ostream& out);

    std::pair<double, double> FreedmanDiaconisBinSize() const;

    size_t size() const;

    bool empty() const;

    QuantileStats<T> get1() const;
    QuantileStats<T> get2() const;

    QuantileStats<T> get(const std::string name) const;

    std::tuple<double, double, double, double> getQuantileEllipse(double const ignore) const;
    void getQuantileEllipse(Ellipses& ellipse, double const ignore) const;

private:

    mutable std::vector<std::pair<T, T>> values;

    QuantileStats<T> quantiles_1;
    QuantileStats<T> quantiles_2;
};

template<class T>
class StatsN {
private:
    std::mutex push_mutex;

    std::vector<std::string> names;

    mutable std::vector<std::vector<T>> values;
public:
    StatsN(std::vector<std::string> const _names);
    StatsN(StatsN<T> const& other);
    StatsN<T> & operator =(StatsN<T> const& other);

    size_t dim() const;

    template<class U>
    bool push(std::vector<U> const& val);

    template<class U>
    bool push_unsafe(std::vector<U> const& val);

    Stats2D<T> getStats2D(size_t ii, size_t jj) const;

    void plotAll(std::string const prefix, HistConfig const conf) const;

    size_t size() const;

    void reserve(const size_t size);
};

template<class T>
class Image2D;

template<class T>
class Image1D {
protected:
    std::vector<T> pos;
    std::vector<T> neg;

    double const width;
private:
    Image2D<T> * parent = nullptr;
public:
    Image1D(double const _width, Image2D<T> * _parent = nullptr);

    T& operator[](double index);

    double min_val = std::numeric_limits<double>::max();
    double max_val = -std::numeric_limits<double>::max();

    void plot(std::string const& prefix, HistConfig const& conf);

    void data2file(std::ostream &out, const HistConfig &conf);
};

template<class T>
class Image2D {
    double const width1;
    double const width2;

    std::vector<Image1D<T> > pos;
    std::vector<Image1D<T> > neg;

public:
    Image2D(double const _width1, const double _width2);

    Image1D<T>& operator[](double const index);

    double min_x = std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();

    double max_x = -std::numeric_limits<double>::max();
    double max_y = -std::numeric_limits<double>::max();

    void plot(std::string const& prefix, HistConfig const& conf);

    void push_unsafe(double const x, double const y, std::vector<double> const& values);

    void data2file(std::ostream& out, HistConfig const& conf);
};

} // namespace runningstats

#endif // RUNNINGSTATS_H
