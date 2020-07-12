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


}
