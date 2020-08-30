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

TEST(RunningStats, copy) {
    RunningStats src, copy, copy2;
    src.push_unsafe(1);
    copy = src;
    EXPECT_EQ(1, copy.getCount());
    EXPECT_EQ(1, src.getCount());
    copy.push_unsafe(1);
    EXPECT_EQ(2, copy.getCount());
    EXPECT_EQ(1, src.getCount());
    copy2 = copy;
    EXPECT_EQ(2, copy.getCount());
    EXPECT_EQ(2, copy2.getCount());
    EXPECT_EQ(1, src.getCount());
    copy2.push_unsafe(1);
    EXPECT_EQ(2, copy.getCount());
    EXPECT_EQ(3, copy2.getCount());
    EXPECT_EQ(1, src.getCount());

    EXPECT_NEAR(1, src.getMean(), 1e-14);
    EXPECT_NEAR(1, copy.getMean(), 1e-14);
    EXPECT_NEAR(1, copy2.getMean(), 1e-14);

    EXPECT_NEAR(0, src.getVar(), 1e-14);
    EXPECT_NEAR(0, copy.getVar(), 1e-14);
    EXPECT_NEAR(0, copy2.getVar(), 1e-14);

    EXPECT_NEAR(0, src.getLogMean(), 1e-14);
    EXPECT_NEAR(0, copy.getLogMean(), 1e-14);
    EXPECT_NEAR(0, copy2.getLogMean(), 1e-14);

    EXPECT_NEAR(0, src.getLogVar(), 1e-14);
    EXPECT_NEAR(0, copy.getLogVar(), 1e-14);
    EXPECT_NEAR(0, copy2.getLogVar(), 1e-14);
}

TEST(StatsN, copy) {
    StatsN<float> src({"a","b"}), copy({"a","b"});
    src.push_unsafe<float>({1,2});
    copy = src;
    EXPECT_EQ(1, copy.size());
    EXPECT_EQ(1, src.size());
    copy.push_unsafe<float>({1,2});
    EXPECT_EQ(2, copy.size());
    EXPECT_EQ(1, src.size());
    StatsN<float> copy2(copy);
    EXPECT_EQ(2, copy.size());
    EXPECT_EQ(2, copy2.size());
    EXPECT_EQ(1, src.size());
    copy2.push_unsafe<float>({1,2});
    EXPECT_EQ(2, copy.size());
    EXPECT_EQ(3, copy2.size());
    EXPECT_EQ(1, src.size());

}
