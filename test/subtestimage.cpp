#include <iostream>
#include <random>

#include <gtest/gtest.h>

#include "runningstats/runningstats.h"

using namespace runningstats;

namespace {
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
}

TEST(Image1D, double) {
    {
        Image1D<double> test(1);
        for (int ii = 0; ii < 10; ++ii) {
            test[ii]++;
        }
        for (int ii = 0; ii < 10; ++ii) {
            EXPECT_NEAR(1, test[ii], 1e-15);
        }
    }
    {
        Image1D<double> test(.25);
        for (double ii = .5; ii < 10; ii += .5) {
            test[ii] = 15.3;
        }
        for (double ii = .5; ii < 10; ii += .5) {
            EXPECT_NEAR(15.3, test[ii], 1e-14);
        }
        for (double ii = .25; ii < 10; ii += .5) {
            EXPECT_NEAR(0, test[ii], 1e-14);
        }
        for (double ii = -.5; ii >= -10; ii -= .5) {
            test[ii] = -15.3;
        }
        for (double ii = -10; ii < 10; ii += .5) {
            EXPECT_NEAR(sgn(ii) * 15.3, test[ii], 1e-14);
        }
        for (double ii = -10.25; ii < 11; ii += .5) {
            EXPECT_NEAR(0, test[ii], 1e-14);
        }
    }
}

TEST(Image1D, RunningStats) {
    {
        Image1D<RunningStats> test(1);
        for (int ii = 0; ii < 10; ++ii) {
            test[ii].push_unsafe(0);
            test[ii].push_unsafe(-1);
            test[ii].push_unsafe(1);
        }
        for (int ii = 0; ii < 10; ++ii) {
            EXPECT_NEAR(0, test[ii].getMean(), 1e-15);
            EXPECT_NEAR(1, test[ii].getStddev(), 1e-15);
        }
    }
    {
        Image1D<RunningStats> test(.25);
        for (double ii = 0; ii < 10; ii += .5) {
            test[ii].push_unsafe(0);
            test[ii].push_unsafe(-1);
            test[ii].push_unsafe(1);
        }
        for (double ii = 0; ii < 10; ii += .5) {
            EXPECT_NEAR(0, test[ii].getMean(), 1e-14);
            EXPECT_NEAR(1, test[ii].getStddev(), 1e-15);
        }
        for (double ii = .25; ii < 10; ii += .5) {
            EXPECT_NEAR(0, test[ii].getMean(), 1e-14);
            EXPECT_NEAR(0, test[ii].getStddev(), 1e-15);
        }
        for (double ii = -.5; ii >= -10; ii -= .5) {
            test[ii].push_unsafe(0);
            test[ii].push_unsafe(-1);
            test[ii].push_unsafe(1);
        }
        for (double ii = -10; ii < 10; ii += .5) {
            EXPECT_NEAR(0, test[ii].getMean(), 1e-14);
            EXPECT_NEAR(1, test[ii].getStddev(), 1e-15);
        }
        for (double ii = -10.25; ii < 11; ii += .5) {
            EXPECT_NEAR(0, test[ii].getMean(), 1e-14);
            EXPECT_NEAR(0, test[ii].getStddev(), 1e-15);
        }
    }
}
