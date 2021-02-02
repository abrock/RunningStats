#include <iostream>
#include <random>

#include <gtest/gtest.h>

#include "runningstats/runningstats.h"

using namespace runningstats;

#include "randutils.hpp"

TEST(ThresholdErrorMean, simple) {
    ThresholdErrorMean<float> stats;
    for (size_t ii = 0; ii < 10000; ++ii) {
        stats.push_unsafe(double(ii) * sqrt(ii), std::sqrt(ii));
    }
    stats.plot("error-threshold-nonrandom", HistConfig().setXLabel("Confidence measure").setMaxPlotPts(1000));
}

TEST(ThresholdErrorMean, random) {
    ThresholdErrorMean<float> stats;
    randutils::mt19937_rng rng;
    for (size_t ii = 0; ii < 10000; ++ii) {
        float const true_error = rng.variate<float, std::normal_distribution>(0,1);
        float const error_measure = true_error + rng.variate<float, std::normal_distribution>(0,1);
        stats.push_unsafe(error_measure, true_error);
    }
    stats.plot("error-threshold-random", HistConfig().setXLabel("Confidence measure").setMaxPlotPts(1000));
}

#include "randutils.hpp"

TEST(ThresholdErrorMean, save_load) {
    ThresholdErrorMean<float> stats;
    for (size_t ii = 0; ii < 100; ++ii) {
        stats.push_unsafe(double(ii) * sqrt(ii), std::sqrt(ii));
    }
    std::string const filename = "threshold-error-mean-test1.bin";
    stats.save(filename);
    ThresholdErrorMean<float> stats2;
    stats2.load(filename);
    ASSERT_EQ(stats.size(), stats2.size());
    ASSERT_EQ(stats.size(), 100);
    for (size_t ii = 0; ii < stats.size(); ++ii) {
        ASSERT_EQ(stats.getData()[ii].first,  stats2.getData()[ii].first);
        ASSERT_EQ(stats.getData()[ii].second, stats2.getData()[ii].second);
    }

}

TEST(ThresholdErrorMean, save_load_random) {
    ThresholdErrorMean<float> stats;
    randutils::mt19937_rng gen;
    size_t const num_test = 1'000'000;
    for (size_t ii = 0; ii < num_test; ++ii) {
        stats.push_unsafe(gen.variate<float, std::normal_distribution>(0,1),
                          gen.variate<float, std::normal_distribution>(0,1));
    }
    std::cout << "Filled stats" << std::endl;
    ASSERT_EQ(stats.size(), num_test);
    std::string const filename = "threshold-error-mean-test1.bin";

#pragma omp parallel sections
    {
#pragma omp section
        stats.save(filename);
#pragma omp section
        stats.plot("threshold-error-mean-test1-plots", HistConfig()
                   .setMaxPlotPts(200));
    }
    std::cout << "Saved stats" << std::endl;
    ThresholdErrorMean<float> stats2;
    stats2.load(filename);
    std::cout << "Loaded stats" << std::endl;
    ASSERT_EQ(stats.size(), stats2.size());
    ASSERT_EQ(stats2.size(), num_test);
    std::vector<std::pair<float, float> > d1 = stats.getData();
    std::vector<std::pair<float, float> > d2 = stats2.getData();

    for (size_t ii = 0; ii < num_test; ++ii) {
        ASSERT_EQ(d1[ii].first,  d2[ii].first);
        ASSERT_EQ(d1[ii].second, d2[ii].second);
    }
    std::cout << "Compared stats" << std::endl;
}
