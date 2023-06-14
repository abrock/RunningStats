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
    disable_thousands_separator(cmd);
    std::string data_file = prefix + ".data";
    Misc::make_target_dir(data_file);
    double const size_x = limits_x.getMax() - limits_x.getMin();
    double const size_y = limits_y.getMax() - limits_y.getMin();
    double const margin = 0.04;
    cmd << "#!/usr/bin/gnuplot \n";
    cmd << "set term png;\n"
        << "set output '" << prefix + ".png'; \n"
        << conf.toString() << ";\n"
        << "set xrange [" << limits_x.getMin() - margin * size_x << ":" << limits_x.getMax() + margin * size_x << "];\n"
        << "set yrange [" << limits_y.getMin() - margin * size_y << ":" << limits_y.getMax() + margin * size_y << "];\n"
           ;

    cmd << "plot " << gpl.file(data, data_file) << " with ellipses notitle; \n";
    cmd << conf.generateFormatCommands(prefix);
    //cmd << "set term tikz; \n"
    //<< "set output '" << prefix << ".tex'; \n"
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

    std::vector<HistConfig::ExtractConf> extractors = conf.extractors;
    if (extractors.empty()) {
        extractors.push_back({conf.extract, conf.extractParam, "", ""});
    }

    std::string const error_file1 = prefix + "-errors-over-threshold";
    std::string const percentage_file1 = prefix + "-percentages-over-threshold";

    gnuplotio::Gnuplot plt1;
    Misc::make_target_dir(prefix);
    std::stringstream cmd1;
    disable_thousands_separator(cmd1);
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
    disable_thousands_separator(cmd2);
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
    for (HistConfig::ExtractConf const& extractor : extractors) {
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
                double const error = error_stats.getStat(extractor.ex, extractor.value);
                double const oracle_error = oracle_stats.getStat(extractor.ex, extractor.value);
                double const threshold = it.first;
                double const percentage = 100.0*double(ii) / double(data.size());
                error_over_threshold.push_back({threshold, error});
                percentage_over_threshold.push_back({threshold, percentage});
                error_over_percentage.push_back({percentage, error});
                error_over_percentage_oracle.push_back({percentage, oracle_error});
                threshold_over_percentage.push_back({percentage, threshold});
            }
        }
        double const error = error_stats.getStat(extractor.ex, extractor.value);
        double const threshold = data.rbegin()->first;
        double const percentage = 100.0;
        error_over_threshold.push_back({threshold, error});
        percentage_over_threshold.push_back({threshold, percentage});
        error_over_percentage.push_back({percentage, error});
        threshold_over_percentage.push_back({percentage, threshold});
        std::string const extractName = extractor.getName();
        cmd1 << plt1.file(error_over_threshold, error_file1 + extractName + ".data")
             << " u 1:2 with l title 'error " << escape(extractName) << "', ";
        cmd2 << plt2.file(error_over_percentage, error_file2 + extractName + ".data")
             << " u 1:2 with l title 'error " << escape(extractName) << "',";
        cmd2 << plt2.file(error_over_percentage_oracle, error_file2 + extractName + "-oracle.data")
             << " u 1:2 with l title 'oracle " << escape(extractName) << "',";

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

std::string escape(const std::string &str) {
    std::string res;
    for (const char c : str) {
        if ('\'' == c) {
            res += '\'';
        }
        res += c;
    }
    return res;
}

std::string tostring(int64_t n) {
    std::string result = std::to_string(n);
    if (result.size() < 4) {
        return result;
    }
    std::string _res;
    _res += result[0];
    for (size_t ii = 1; ii < result.size(); ++ii) {
        if (ii+1 < result.size()) {
            if ((result.size() - ii) % 3 == 0) {
                _res += '\'';
            }
        }
        _res += result[ii];
    }
    return _res;
}

void disable_thousands_separator(std::ostream &out) {
    std::locale loc("");
    // imbue loc and add your own facet:
    out.imbue(std::locale(loc, new no_separator()));
}

} // namespace runningstats
