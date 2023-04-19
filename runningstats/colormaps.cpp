#include "runningstats.h"

namespace runningstats {

ColorMaps::ColorMaps() {
#include "colormaps/blue_red_lab_ceres.inc"
#include "colormaps/blue_red_lab_ceres_clipped.inc"
#include "colormaps/blue_red_lab_ceres2.inc"
#include "colormaps/viridis.inc"
#include "colormaps/viridis-clipped.inc"
}

std::string ColorMaps::blue_red() {
    return "blue-red";
}

std::string ColorMaps::blue_red_clipped() {
    return "blue-red-clipped";
}

std::string ColorMaps::blue_red_2() {
    return "blue-red-2";
}

std::string ColorMaps::viridis() {
    return "viridis";
}

std::string ColorMaps::viridis_clipped() {
    return "viridis-clipped";
}

}
