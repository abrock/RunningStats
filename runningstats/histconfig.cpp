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

void LineSegment::addPt(const float x, const float y) {
    data.push_back({x,y});
}

void Line::addSegment(LineSegment const& seg) {
    data.push_back(seg);
}

std::vector<std::vector<std::pair<float, float> > > Line::getData() const {
    std::vector<std::vector<std::pair<float, float> > > result;
    for (LineSegment const& line : data) {
        result.push_back(line.data);
    }
    return result;
}

void LineSegment::write(std::ostream& out) const {
    out.precision(15);
    disable_thousands_separator(out);
    for (const std::pair<float, float> &pair : data) {
        out << pair.first << "\t" << pair.second << std::endl;
    }
    out << std::endl;
}

void Line::writeToFile(std::string const& filename) const {
    Misc::make_target_dir(filename);
    std::ofstream out(filename, std::ostream::out);
    for (LineSegment const& seg : data) {
        seg.write(out);
    }
}

std::string Line::getFilename(const std::string &prefix, const std::string & suffix) const {
    return prefix + "-contour-" + suffix + ".data";
}

bool Line::empty() const {
    return data.empty();
}

bool LineSegment::empty() const {
    return data.empty();
}

void Line::clear() {
    data.clear();
}

void LineSegment::clear() {
    data.clear();
}

std::string Contour::generateContour(const std::string &data_file) const {
    Misc::make_target_dir(data_file);
    return "\nset contour\n"
    "unset surface\n"
    "set cntrparam levels discrete " + std::to_string(value) + "\n"
    "set view map\n"
    "unset clabel\n"
    "set table '" + contoursFilename(data_file) + "'\n"
    "splot '" + data_file + "' u 1:2:3 not\n"
    "unset table\n"
    "unset contour\n";
}

std::string Contour::plotContour(const std::string &data_file) const {
    Misc::make_target_dir(data_file);
    return
            ", '" + contoursFilename(data_file) + "' w l " +
            (color.empty() ? "" : "lc rgb '#" + color + "' ") +
            (title.empty() ? "notitle " : "title '" + title + "' ");
}

std::string Contour::contoursFilename(const std::string &prefix) const {
    return prefix + "-contour-" + std::to_string(value);
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

HistConfig &HistConfig::addContour(const double value, const std::string _title, const std::string color) {
    contours[value] = {value, _title, color};
    return *this;
}

HistConfig &HistConfig::addLine(Line const& l) {
    lines.push_back(l);
    return *this;
}


std::string HistConfig::generateContours(const std::string &data_file) const {
    std::string result;
    for (auto const& it : contours) {
        result += it.second.generateContour(data_file);
    }
    return result;
}

std::string HistConfig::createTics(const std::vector<double> &positions, const std::vector<std::string> &names) {
    if (positions.size() != names.size()) {
        throw std::runtime_error("Number of elements in position and name vectors must be the same but is different: " +
                                 std::to_string(positions.size()) + " vs. " + std::to_string(names.size()));
    }
    if (positions.empty()) {
        return "";
    }
    int mod = 1;
    if (positions.size() > maxLabeledTics) {
        mod = std::ceil(double(positions.size()) / maxLabeledTics);
    }
    std::stringstream result;
    for (size_t ii = 0; ii < positions.size(); ++ii) {
        result << '"' << (ii % mod == 0 ? names[ii] : "") << '"' << " " << positions[ii] << ", ";
    }
    return result.str();
}

std::string HistConfig::plotContours(const std::string &data_file) const {
    std::string result;
    for (auto const& it : contours) {
        result += it.second.plotContour(data_file);
    }
    return result;
}

std::string HistConfig::plotLines(const std::string &prefix) const {
    std::stringstream result;
    disable_thousands_separator(result);
    for (size_t ii = 0; ii < lines.size(); ++ii) {
        Line const& l = lines[ii];
        const std::string _title = l.title.empty() ? " notitle " : " title '" + escape(l.title) + "' ";
        const std::string color = l.color.empty() ? " " : " lc rgb '" + escape(l.color) + "' ";
        const std::string fn = l.getFilename(prefix, "-line-" + std::to_string(ii) + l.title);
        l.writeToFile(fn);
        result << ", '" << escape(fn) << "'"
            << " u 1:2 w l " << color << _title;
    }
    return result.str();
}

HistConfig& HistConfig::setMaxPlotPts(int64_t val) {
    max_plot_pts = val;
    return *this;
}

HistConfig &HistConfig::setIgnoreAmount(const double val) {
    ignore_amount = std::min(1.0, std::max(0.0, val));
    return *this;
}

std::string HistConfig::colorMapCmd() const {
    std::stringstream cmd;
    if (!colormap.empty()) {
        std::string const pal_name = colormap + ".pal";
        cmd << "load '" << pal_name << "';\n";
        static std::mutex lock;
        std::lock_guard guard(lock);
        if (!fs::exists(pal_name)) {
            std::ofstream colormap_out(pal_name);
            colormap_out << ColorMaps().getColorMap(colormap) << std::endl;
        }
    }
    return cmd.str();
}

HistConfig& HistConfig::setBG(const std::string &color) {
    bg_color = color;
    return *this;
}

void HistConfig::getLinesRect(double &_min_x, double &_min_y, double &_max_x, double &_max_y) const {
    if (lines.empty()) {
        return;
    }
    bool found_first = false;
    for (Line const& l: lines) {
        for (LineSegment const& s : l.data) {
            for (std::pair<float, float> const& pt : s.data) {
                if (!found_first) {
                    _min_x = _max_x = pt.first;
                    _min_y = _max_y = pt.second;
                    found_first = true;
                }
                _min_x = std::min<double>(_min_x, pt.first);
                _max_x = std::max<double>(_max_x, pt.first);

                _min_y = std::min<double>(_min_y, pt.second);
                _max_y = std::max<double>(_max_y, pt.second);
            }
        }
    }
}

HistConfig &HistConfig::setColorMap(const std::string &name) {
    ColorMaps().getColorMap(name);
    colormap = name;
    return *this;
}

std::string HistConfig::toString() const {
    std::stringstream out;
    disable_thousands_separator(out);
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
    if (!xTics.empty()) {
        out << "set xtics (" << xTics << ");\n";
    }
    if (!yTics.empty()) {
        out << "set ytics (" << yTics << ");\n";
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
    case Extract::Count: return "Count";
    }
    return "Unknown";
}

std::string HistConfig::ExtractConf::getName() const {
    switch (ex) {
    case Extract::Mean : return extractName(ex);
    case Extract::Median : return extractName(ex);
    case Extract::Stddev : return extractName(ex);
    case Extract::Variance : return extractName(ex);
    case Extract::Count: return extractName(ex);
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
    if (s == "count") return Extract::Count;
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
    disable_thousands_separator(cmd);
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

HistConfig &HistConfig::extractCount() {
    extract = Extract::Count;
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

HistConfig &HistConfig::setMinMaxCB(double _min_cb, double _max_cb) {
    this->min_cb = _min_cb;
    this->max_cb = _max_cb;
    return *this;
}

HistConfig &HistConfig::setMinCB(double _min_cb) {
    this->min_cb = _min_cb;
    return *this;
}

HistConfig &HistConfig::setMaxCB(double _max_cb) {
    this->max_cb = _max_cb;
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


HistConfig& HistConfig::setXtics(const std::vector<double> &positions, const std::vector<std::string> &names) {
    xTics = createTics(positions, names);
    return *this;
}

HistConfig& HistConfig::setYtics(const std::vector<double> &positions, const std::vector<std::string> &names) {
    yTics = createTics(positions, names);
    return *this;
}

std::string HistConfig::cbRangeToString() const {
    return "set cbrange " + Misc::range2string(min_cb, max_cb) + ";\n";
}

template<class T>
HistConfig& HistConfig::setSymmetricCBRange(QuantileStats<T> const& stats, double quantile) {
    quantile = std::min(1.0, std::max(0.0, quantile));
    double const range = std::max(stats.getQuantile(quantile), stats.getQuantile(1.0 - quantile));
    setMinMaxCB(-range, range);
    return *this;
}

template HistConfig& HistConfig::setSymmetricCBRange(QuantileStats<float> const&, double);
template HistConfig& HistConfig::setSymmetricCBRange(QuantileStats<double> const&, double);

template<class T>
HistConfig& HistConfig::setSymmetricCBRange(Image2D<T> const& stats, double quantile) {
    setSymmetricCBRange(stats.template merged<float>(), quantile);
    return *this;
}

template HistConfig& HistConfig::setSymmetricCBRange(Image2D<float> const&, double);
template HistConfig& HistConfig::setSymmetricCBRange(Image2D<QuantileStats<float> > const&, double);

} // namespace runningstats
