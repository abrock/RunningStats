#include <iostream>
#include <random>

#include <gtest/gtest.h>

#include "runningstats/runningstats.h"

using namespace runningstats;

static std::random_device dev;
static std::mt19937 engine(dev());

TEST(Plot2D, Histogram2D) {
    std::normal_distribution<double> dist;

    runningstats::Histogram2D hist(0.1, 0.2);

    for (size_t ii = 0; ii < 1000*1000; ++ii) {
        hist.push_unsafe(.2*dist(engine), dist(engine));
    }

    hist.plotHist("normal-2d", false);

}

TEST(Plot2D, Stats2D) {
    std::normal_distribution<double> dist;

    runningstats::Stats2D<float> stats;

    for (size_t ii = 0; ii < 1000*1000; ++ii) {
        stats.push_unsafe(.1*dist(engine), dist(engine));
    }

    runningstats::Histogram2D hist = stats.getHistogram2D(stats.FreedmanDiaconisBinSize());
    hist.plotHist("stats-normal-2d", false);

}

TEST(Plot2D, Stats2Dfixed) {
    std::normal_distribution<double> dist;

    runningstats::Stats2D<float> stats;

    for (size_t ii = 0; ii < 1000*1000; ++ii) {
        stats.push_unsafe(.1*dist(engine), dist(engine));
    }

    runningstats::Histogram2Dfixed hist = stats.getHistogram2Dfixed(stats.FreedmanDiaconisBinSize());
    hist.plotHist("stats-fixed-normal-2d", false);

}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    std::cout << "RUN_ALL_TESTS return value: " << RUN_ALL_TESTS() << std::endl;
}
