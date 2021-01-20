#include <iostream>
#include <random>

#include <gtest/gtest.h>

#include "runningstats/runningstats.h"

using namespace runningstats;

TEST(WeightedRunningStats, mean) {
    for (int limit = 0; limit < 1'000; ++limit) {
        WeightedRunningStats s;
        for (int ii = -limit; ii <= limit; ++ii) {
            s.push_unsafe(ii, std::abs(ii));
        }
        EXPECT_NEAR(s.getMean(), 0, 1e-12);
        for (int ii = 0; ii <= limit; ++ii) {
            s.push_unsafe(ii, 0);
        }
        EXPECT_NEAR(s.getMean(), 0, 1e-12);
    }
}

TEST(WeightedRunningStats, var) {
    for (int limit = 0; limit < 1'000; ++limit) {
        WeightedRunningStats s;
        for (int ii = 0; ii <= limit; ++ii) {
            s.push_unsafe(1, 1);
            s.push_unsafe(-1, 1);
        }
        EXPECT_NEAR(s.getMean(), 0, 1e-12);
        EXPECT_NEAR(s.getVar(), 1, 1e-12);
    }
    for (int limit = 1; limit < 1'000; ++limit) {
        WeightedRunningStats s;
        for (int ii = 0; ii <= limit; ++ii) {
            s.push_unsafe(1, ii);
            s.push_unsafe(-1, ii);
        }
        EXPECT_NEAR(s.getMean(), 0, 1e-12);
        EXPECT_NEAR(s.getVar(), 1, 1e-12);
    }

    for (int low_val = -10; low_val <= 10; ++low_val) {
        for (int high_val = -10; high_val <= 10; ++high_val) {
            WeightedRunningStats s;
            for (int ii = 0; ii <= 20; ++ii) {
                s.push_unsafe(low_val, ii);
                s.push_unsafe(high_val, ii);
            }
            EXPECT_NEAR(s.getMean(), double(low_val + high_val)/2, 1e-12);
            double const diff = std::abs(double(low_val - high_val))/2;
            EXPECT_NEAR(s.getVar(), diff*diff, 1e-12);
        }
    }

    for (int mean = -3; mean <= 3; ++mean) {
        for (int max_offset = 1; max_offset < 10; ++max_offset) {
            WeightedRunningStats s;
            double squaresum = 0;
            for (int offset = 0; offset < max_offset; ++offset) {
                squaresum += offset*offset;
                for (double weight = 0.25; weight < 20; weight *= 1.25) {
                    s.push_unsafe(mean + offset, weight);
                    s.push_unsafe(mean - offset, weight);
                }
            }
            EXPECT_NEAR(s.getVar(), squaresum/max_offset, 1e-12);
        }
    }
}
