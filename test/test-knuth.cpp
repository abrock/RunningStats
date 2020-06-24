#include <iostream>
#include <random>
#include <fstream>

#include <gtest/gtest.h>

#include "runningstats/runningstats.h"

using namespace runningstats;

static std::random_device dev;
static std::mt19937 engine(dev());

TEST(Knuth, brute_force) {
    runningstats::QuantileStats<float> data;
    std::normal_distribution<double> dist;
    //std::uniform_real_distribution<double> dist;
    for (size_t ii = 0; ii < 1*1000; ++ii) {
        data.push_unsafe(dist(engine));
    }
    double const range = data.getMax() - data.getMin();
    std::ofstream data_out("normal-likelihood.data");
    double best_val = -std::numeric_limits<double>::max();
    size_t best_M = 0;
    size_t max_M = 4*range/data.FreedmanDiaconisBinSize();
    for (size_t M = 1; M < max_M; ++M) {
        double const bin_width = range / M;
        runningstats::Histogram hist = data.getHistogram(bin_width);
        double const prob = hist.getPosteriorProbability();
        data_out << M << "\t" << prob << std::endl;
        if (prob >= best_val) {
            best_M = M;
            best_val = prob;
        }
    }
    std::cout << "Max M: " << max_M << std::endl;
    std::cout << "Best M: " << best_M << std::endl
              << "Best bin width: " << range / best_M << std::endl;
    data.plotHist("normal-knuth", range / best_M, false);
    data.plotHist("normal-freedman-diaconis", data.FreedmanDiaconisBinSize(), false);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    std::cout << "RUN_ALL_TESTS return value: " << RUN_ALL_TESTS() << std::endl;
}
