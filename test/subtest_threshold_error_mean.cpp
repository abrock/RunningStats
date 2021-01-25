#include <iostream>
#include <random>

#include <gtest/gtest.h>

#include "runningstats/runningstats.h"

using namespace runningstats;

TEST(ThresholdErrorMean, simple) {
    ThresholdErrorMean<float> stats;
    for (size_t ii = 0; ii < 10000; ++ii) {
        stats.push_unsafe(double(ii) * sqrt(ii), std::sqrt(ii));
    }

    stats.plot("error-threshold", HistConfig().setXLabel("Confidence measure").setMaxPlotPts(1000));



}
