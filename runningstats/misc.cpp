#include "misc.h"

#include <cmath>

#include <filesystem>
namespace fs = std::filesystem;

namespace {
std::error_code ignore_error_code;
}


void Misc::make_target_dir(const std::string &filename) {
    fs::create_directories(fs::path(filename).remove_filename(), ignore_error_code);
}

std::string Misc::range2string(const double min, const double max) {
    std::string const _min = std::isfinite(min) ? std::to_string(min) : "";
    std::string const _max = std::isfinite(max) ? std::to_string(max) : "";
    return "[" + _min + ":" + _max + "]";
}

std::string Misc::replace_empty(const std::string &val, const std::string &replacement) {
    return val.empty() ? replacement : val;
}

