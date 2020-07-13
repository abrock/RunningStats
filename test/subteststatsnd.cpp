#include <iostream>
#include <random>

#include <gtest/gtest.h>

#include "runningstats/runningstats.h"

using namespace runningstats;

static std::random_device dev;
static std::mt19937 engine(dev());

TEST(StatsN, main) {
    StatsN<float> stats(3, {"A", "B", "A+B"});
    std::normal_distribution<double> dist;
    for (size_t ii = 0; ii < 1000*1000; ++ii) {
        double const A = dist(engine);
        double const B = dist(engine);
        stats.push<double>({A,B,A+B});
    }

    stats.plotAll("statsN-plot-all", HistConfig().setLogCB());
}
