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

#include <filesystem>
namespace fs = std::filesystem;

namespace runningstats {

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
void Image1D<T>::setParent(Image2D<T> &p) {
    parent = &p;
}

template<class T>
T Image1D<T>::aggregate(std::function<T (T, T)> func) const {
    if (neg.empty()) {
        return aggregate(pos, func);
    }
    if (pos.empty()) {
        return aggregate(neg, func);
    }
    return func(aggregate(pos, func), aggregate(neg, func));
}

template<class T>
bool Image1D<T>::empty() const {
    return pos.empty() && neg.empty();
}

template<class T>
T Image1D<T>::aggregate(const std::vector<T> &vec, std::function<T (T, T)> func) {
    if (vec.empty()) {
        return T();
    }
    T result = vec.front();
    if (vec.size() == 1) {
        return result;
    }
    for (size_t ii = 1; ii < vec.size(); ++ii) {
        result = func(result, vec[ii]);
    }
    return result;
}

template<class T>
void Image1D<T>::plot(const std::string &prefix, const HistConfig &conf) {
    std::string const data_file = prefix + ".data";
    Misc::make_target_dir(data_file);
    std::ofstream data_out(data_file);
    data2file(data_out, conf);
    std::stringstream cmd;
    cmd << "#!/usr/bin/gnuplot \n";
    cmd << "set term png;\n";
    cmd << "set output '" << prefix << ".png';\n";
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
    cmd << conf.generateFormatCommands(prefix);
    gnuplotio::Gnuplot plt;
    plt << cmd.str();
    std::ofstream cmd_out(prefix + ".gpl");
    cmd_out << cmd.str();
}

template<class T>
Image2D<T>::Image2D(const double _width1, double const _width2) : width1(_width1), width2(_width2) {
}

template<class T>
Image2D<T>::Image2D(const Image2D<T> &other) {
    width1 = other.width1;
    width2 = other.width2;
    pos = other.pos;
    neg = other.neg;
    for (Image1D<T> & it : pos) {
        it.setParent(*this);
    }
    for (Image1D<T> & it : neg) {
        it.setParent(*this);
    }
}

template<class T>
Image1D<T> &Image2D<T>::operator[](const double index) {
    static std::mutex enlarge_mutex;
    std::lock_guard guard(enlarge_mutex);
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
    Misc::make_target_dir(data_file);
    std::ofstream data_out(data_file);
    data2file(data_out, conf);
    std::stringstream cmd;
    cmd << "#!/usr/bin/gnuplot \n";
    cmd << "set term png;\n";
    cmd << "set output '" << prefix << ".png';\n";
    cmd << conf.toString() << "\n";
    cmd << conf.colorMapCmd();
    double _min_x = min_x;
    double _min_y = min_y;
    double _max_x = max_x;
    double _max_y = max_y;
    conf.getLinesRect(_min_x, _min_y, _max_x, _max_y);
    _min_x = std::min(_min_x, min_x) - width1/2;
    _min_y = std::min(_min_y, min_y) - width2/2;
    _max_x = std::max(_max_x, max_x) - width1/2;
    _max_y = std::max(_max_y, max_y) - width2/2;
    if (conf.flip_x) {
        cmd << "set xrange[" << _max_x << ":" << _min_x << "]; \n";
    }
    else {
        cmd << "set xrange[" << _min_x << ":" << _max_x << "]; \n";
    }
    if (conf.flip_y) {
        cmd << "set yrange[" << _max_y << ":" << _min_y << "]; \n";
    }
    else {
        cmd << "set yrange[" << _min_y << ":" << _max_y << "]; \n";
    }
    cmd << "set object rectangle from " << _min_x << "," << _min_y << " to " << _max_x << "," << _max_y <<
           " behind fillcolor rgb '" << escape(conf.bg_color) << "' fillstyle solid noborder;\n";
    cmd << "set xtics out;\n";
    cmd << "set ytics out;\n";
    if (conf.min_cb > -std::numeric_limits<double>::max()) {
        cmd << "set cbrange [" << conf.min_cb << ":" << conf.max_cb << "];\n";
    }
    cmd << conf.generateContours(data_file);
    cmd << "plot '" << data_file << "' u 1:2:3 with image notitle " << conf.plotContours(data_file) << conf.plotLines(data_file);
    cmd << ";\n";
    cmd << conf.generateFormatCommands(prefix);
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

template<class T>
Image2D<T> Image2D<T>::mergeImages(
        Image2D<T> &a,
        Image2D<T> &b,
        std::function<T (T, T)> func) {
    Image2D<T> result(a.width1, a.width2);

    for (double x = std::min(a.min_x, b.min_x); x <= std::max(a.max_x, b.max_x); x += a.width1) {
        for (double y = std::min(a.min_y, b.min_y); y <= std::max(a.max_y, b.max_y); y += a.width2) {
            result[x][y] = func(a[x][y], b[x][y]);
        }
    }

    return result;
}

template<class T>
T Image2D<T>::aggregate(std::function<T (T, T)> func) const {
    std::vector<Image1D<T> const* > non_empty_parts;
    for (Image1D<T> const& t : pos) {
        if (!t.empty()) {
            non_empty_parts.push_back(&t);
        }
    }
    for (Image1D<T> const& t : neg) {
        if (!t.empty()) {
            non_empty_parts.push_back(&t);
        }
    }
    if (non_empty_parts.empty()) {
        return T();
    }
    T result = non_empty_parts.front()->aggregate(func);
    if (non_empty_parts.size() == 1) {
        return result;
    }
    for (size_t ii = 1; ii < non_empty_parts.size(); ++ii) {
        result = func(result, non_empty_parts[ii]->aggregate(func));
    }
    return result;
}

template<class T>
template<class U>
void Image1D<T>::merge(QuantileStats<U> &stats) const {
    for (T const& data : pos) {
        stats.push_unsafe(data);
    }
    for (T const& data : neg) {
        stats.push_unsafe(data);
    }
}

template<class T>
Image2D<T>::Image2D() : width1(1), width2(1) {}

template void Image1D<float>::merge(QuantileStats<float> &stats) const;
template void Image1D<QuantileStats<float> >::merge(QuantileStats<float> &stats) const;


template<class T>
template<class U>
QuantileStats<U> Image2D<T>::merged() const {
    QuantileStats<U> result;
    for (Image1D<T> const& data: pos) {
        data.merge(result);
    }
    for (Image1D<T> const& data: neg) {
        data.merge(result);
    }
    return result;
}

template QuantileStats<float> Image2D<float>::merged() const;
template QuantileStats<float> Image2D<QuantileStats<float> >::merged() const;

template class Image1D<double>;
template class Image1D<RunningStats>;
template class Image1D<QuantileStats<float> >;
template class Image2D<double>;
template class Image2D<float>;
template class Image2D<size_t>;
template class Image2D<RunningStats>;
template class Image2D<QuantileStats<float> >;
template class Image2D<std::vector<QuantileStats<float> > >;

} // namespace runningstats
