#include "runningstats.h"

namespace runningstats {

ColorMaps::ColorMaps() {
#include "colormaps/blue_red_lab_ceres.inc"
#include "colormaps/blue_red_lab_ceres2.inc"
}

std::string ColorMaps::blue_red() {
    return "blue-red";
}

std::string ColorMaps::blue_red_2() {
    return "blue-red-2";
}

}
