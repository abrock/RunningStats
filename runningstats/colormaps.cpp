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
    return std::string("blue-red");
}

std::string ColorMaps::blue_red_clipped() {
    return std::string("blue-red-clipped");
}

std::string ColorMaps::blue_red_2() {
    return std::string("blue-red-2");
}

std::string ColorMaps::viridis() {
    return std::string("viridis");
}

std::string ColorMaps::viridis_clipped() {
    return std::string("viridis-clipped");
}

std::string ColorMaps::getColorMap(const std::string &name) {
    auto const& it = colormaps.find(name);
    if (colormaps.end() != it) {
        return it->second;
    }

    throw std::runtime_error("Colormap " + name + " unknown");
}

}
