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
        out << "set xlabel '" << escape(xLabel) << "';\n";
    }
    if (!yLabel.empty()) {
        out << "set ylabel '" << escape(yLabel) << "';\n";
    }
    if (!title.empty()) {
        out << "set title '" << escape(title) << "';\n";
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

std::string HistConfig::ExtractConf::getName() const {
    switch (ex) {
    case Extract::Mean : return extractName(ex);
    case Extract::Median : return extractName(ex);
    case Extract::Stddev : return extractName(ex);
    case Extract::Variance : return extractName(ex);
    default: break;
    }
    std::string result = extractName(ex) + "-" + std::to_string(value);
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

HistConfig &HistConfig::addExtractor(const HistConfig::Extract _ex, const double _val, const std::string &_color, const std::string &_label) {
    extractors.push_back({_ex, _val, _color, _label});
    return *this;
}

HistConfig& HistConfig::addExtractors(const std::vector<std::pair<std::string, double> > &vec) {
    for (auto const& it : vec) {
        addExtractor(str2extract(it.first), it.second, "", "");
    }
    return *this;
}

HistConfig &HistConfig::addCommand(const std::string &cmd) {
    misc += ";\n" + cmd + ";\n";
    return *this;
}

void HistConfig::addFormat(const std::string &fmt) {
    formats.insert(fmt);
}

std::string HistConfig::generateFormatCommands(std::string const& prefix) const {
    std::stringstream cmd;
    for (std::string const& fmt : formats) {
        if ("svg" == fmt) {
            cmd << "set term svg enhanced background rgb 'white';\n";
            cmd << "set output '" << prefix << ".svg';\n"
                << "replot;\n";
        }
    }
    return cmd.str();
}

HistConfig &HistConfig::addExtractors(const std::vector<std::pair<HistConfig::Extract, double> > & vec) {
    for (std::pair<HistConfig::Extract, double> const& it : vec) {
        addExtractor(it.first, it.second, "", "");
    }
    return *this;
}

HistConfig &HistConfig::setExtractors(const std::vector<std::tuple<HistConfig::Extract, double, std::string> > &vec) {
    extractors.clear();
    for (std::tuple<HistConfig::Extract, double, std::string> const& it : vec) {
        addExtractor(std::get<0>(it), std::get<1>(it), std::get<2>(it), "");
    }
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

HistConfig &HistConfig::setFlipX(const bool val) {
    flip_x = val;
    return *this;
}

HistConfig &HistConfig::setFlipY(const bool val) {
    flip_y = val;
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
    misc += std::string("set arrow from graph 0, first ") + std::to_string(y) + " to graph 1, first " + std::to_string(y) + " nohead lc rgb '" + color + "' front;\n";
    return *this;
}

HistConfig &HistConfig::addVerticalLine(const double x, const std::string color) {
    misc += std::string("set arrow from ") + std::to_string(x) + ", graph 0 to " + std::to_string(x) + ", graph 1 nohead lc rgb '" + color + "' front;\n";
    return *this;
}


HistConfig HistConfig::clone() const {
    return HistConfig(*this);
}


} // namespace runningstats