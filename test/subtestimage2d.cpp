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

TEST(Image2D, double) {
    {
        Image2D<double> test(1, 1);
        for (double ii = -10; ii <= 10; ++ii) {
            test[ii][ii] = 1;
        }
        for (double ii = -10; ii <= 10; ++ii) {
            for (double jj = -10; jj <= 10; ++jj) {
                EXPECT_NEAR(ii == jj, test[ii][jj], 1e-16);
            }
        }
    }
}
