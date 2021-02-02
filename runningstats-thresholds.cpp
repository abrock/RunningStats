#include <iostream>

#include "runningstats/runningstats.h"

namespace rs = runningstats;

int main(int argc, char ** argv) {
    runningstats::ThresholdErrorMean<float> stats;

    if (argc < 2) {
        return EXIT_FAILURE;
    }

    std::string const filename = argv[1];

    stats.load(filename);

    rs::HistConfig conf;
    conf.setMaxPlotPts(200);
    conf.addExtractors({
                           {"Mean", 0.5},
                           {"Quantile", 0.5},
                           {"Quantile", 0.75},
                           {"Quantile", 0.90},
                       });

    stats.plot(filename + "-plot", conf);

    return EXIT_SUCCESS;
}
