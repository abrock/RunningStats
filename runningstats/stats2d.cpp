#undef NDEBUG
#include <cassert>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

#include "runningstats.h"

#include "gnuplot-iostream.h"

namespace runningstats {

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
void Stats2D<T>::plotDots(const std::string prefix, const HistConfig &conf) const {
    std::stringstream cmd;
    disable_thousands_separator(cmd);
    gnuplotio::Gnuplot plt;
    std::string const data_fn = plt.file1d(values, prefix + ".data");
    cmd << "#!/usr/bin/gnuplot \n";
    cmd << "set term pngcairo;\n";
    cmd << "set output '" << prefix << ".png';\n";
    cmd << conf.toString() << ";\n";
    if (conf.fit_line) {
        cmd << "f(x) = a*x+b;\n";
        cmd << "fit f(x) " << data_fn << " via a,b;\n";
        cmd << "plot " << data_fn << " title 'data', f(x) w l title 'fit';\n";
    }
    else {
        cmd << "plot " << data_fn << " notitle;\n";
    }
    plt << cmd.str();
    std::ofstream cmd_out(prefix + ".gpl");
    cmd_out << cmd.str();
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
size_t Stats2D<T>::getCount() const {
    return size();
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

} // namespace runningstats
