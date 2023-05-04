#include "misc.h"

#include <filesystem>
namespace fs = std::filesystem;

namespace {
std::error_code ignore_error_code;
}


void Misc::make_target_dir(const std::string &filename) {
    fs::create_directories(fs::path(filename).remove_filename(), ignore_error_code);
}
