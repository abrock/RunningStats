#include <iostream>
#include <random>

#include <gtest/gtest.h>

#include "runningstats/runningstats.h"

using namespace runningstats;

#include "randutils.hpp"

TEST(QuantileStats2, save_load) {
    randutils::mt19937_rng rng;
    size_t const num = 10'000;
    std::string const filename = "QuantileStats2-save-load.bin";

    QuantileStats<float> src, dst;

    for (size_t ii = 0; ii < num; ++ii) {
        src.push_unsafe(rng.variate<float, std::normal_distribution>(0,1));
    }

    src.save(filename);
    dst.load(filename);

    std::vector<float> data1 = src.getData();
    std::vector<float> data2 = dst.getData();
    ASSERT_EQ(data1.size(), data2.size());
    ASSERT_EQ(data1.size(), num);

    ASSERT_EQ(src.getMax(), dst.getMax());
    ASSERT_EQ(src.getMin(), dst.getMin());
    ASSERT_EQ(src.getMean(), dst.getMean());
    ASSERT_EQ(src.getTrimmedMean(.5), dst.getTrimmedMean(.5));
    ASSERT_EQ(src.getQuantile(.5), dst.getQuantile(.5));
    for (double ii = 0; ii <= 1; ii += 1.0/64) {
        ASSERT_EQ(src.getQuantile(ii), dst.getQuantile(ii));
    }

    for (size_t ii = 0; ii < num; ++ii) {
        ASSERT_EQ(data1[ii], data2[ii]);
    }
    src.sort();
    data1 = src.getData();
    float prev = data1.front();
    for (auto const it : data1) {
        ASSERT_LE(prev, it);
        prev = it;
    }
    ASSERT_EQ(src.getMax(), dst.getMax());
    ASSERT_EQ(src.getMin(), dst.getMin());
    ASSERT_EQ(src.getMean(), dst.getMean());
    ASSERT_EQ(src.getTrimmedMean(.5), dst.getTrimmedMean(.5));
    ASSERT_EQ(src.getQuantile(.5), dst.getQuantile(.5));

    src.plotCDF("10k-normal-dist-src");
    dst.plotCDF("10k-normal-dist-dst");
}


