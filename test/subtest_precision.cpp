#include <iostream>
#include <random>

#include <gtest/gtest.h>

#include "runningstats/runningstats.h"

using namespace runningstats;

bool sortbyfirst(std::pair<float, float> const& a, std::pair<float, float> const& b) {
    return a.first < b.first;
}

TEST(Precision, ordering) {
    int const size = 10'000;
    std::vector<std::pair<float, float> > data;
    for (int ii = 0; ii < size; ++ii) {
        data.push_back({-ii, ii});
    }
    ASSERT_EQ(data.size(), size);
    std::sort(data.begin(), data.end(), sortbyfirst);
    for (int ii = 0; ii < size; ++ii) {
        ASSERT_EQ(data[ii].first, -(size - ii - 1));
        ASSERT_EQ(data[ii].second, (size - ii - 1));
    }
}


TEST(Precision, RunningStats) {
    int64_t max_n = 500'000'000;
    RunningStats rs;
    for (int64_t n = 0; n < max_n; ++n) {
        double const val = -pow(1.001, -n);
        rs.push_unsafe(val);
    }
    for (int64_t n = 0; n < max_n; ++n) {
        double const val = pow(1.001, -(max_n-n-1));
        rs.push_unsafe(val);
    }
    std::cout << rs.print() << std::endl;
    // -1.065329840014954e-18 +- 1.000749844327757e-03, 1000000000 Samples, range: [-1.000000000000000e+00, 1.000000000000000e+00
}
