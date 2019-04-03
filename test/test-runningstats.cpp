#include <iostream>

#include <gtest/gtest.h>

#include "runningstats/runningstats.h"

::testing::AssertionResult RelativeNear(const double a, const double b, double delta) {
    double const diff = std::abs(a-b);
    double const relative_diff = 2*diff/(std::abs(a) + std::abs(b));
    if (relative_diff < delta)
        return ::testing::AssertionSuccess();
    else
        return ::testing::AssertionFailure() << "The absolute difference between "
                                             << a << " and " << b << " is " << diff
                                             << ", the relative difference is " << relative_diff
                                             << " which exceeds " << delta;
}

TEST(BinaryStats, Everything) {
    BinaryStats a;
    for (size_t ii = 0; ii < 90; ++ii) {
        a.push(true);
    }
    EXPECT_NEAR(a.get(), 1.0, 1e-15);
    EXPECT_NEAR(a.getPercent(), 100.0, 1e-15);
    EXPECT_EQ(a.getTotalCount(), 90);
    EXPECT_EQ(a.getTrueCount(), 90);
    EXPECT_EQ(a.getFalseCount(), 0);
    for (size_t ii = 0; ii < 10; ++ii) {
        a.push(false);
    }
    EXPECT_NEAR(a.get(), .9, 1e-15);
    EXPECT_NEAR(a.getPercent(), 90.0, 1e-15);
    EXPECT_EQ(a.getTotalCount(), 100);
    EXPECT_EQ(a.getTrueCount(), 90);
    EXPECT_EQ(a.getFalseCount(), 10);
    for (size_t ii = 0; ii < 80; ++ii) {
        a.push(false);
    }
    EXPECT_NEAR(a.get(), .5, 1e-15);
    EXPECT_NEAR(a.getPercent(), 50.0, 1e-15);
    EXPECT_EQ(a.getTotalCount(), 180);
    EXPECT_EQ(a.getTrueCount(), 90);
    EXPECT_EQ(a.getFalseCount(), 90);

    a.pushFalse(10);
    a.pushTrue(10);
    EXPECT_NEAR(a.get(), .5, 1e-15);
    EXPECT_NEAR(a.getPercent(), 50.0, 1e-15);
    EXPECT_EQ(a.getTotalCount(), 200);
    EXPECT_EQ(a.getTrueCount(), 100);
    EXPECT_EQ(a.getFalseCount(), 100);

    BinaryStats b;
    b.pushTrue(800);

    EXPECT_NEAR(b.get(), 1.0, 1e-15);
    EXPECT_NEAR(b.getPercent(), 100.0, 1e-15);
    EXPECT_EQ(b.getTotalCount(), 800);
    EXPECT_EQ(b.getTrueCount(), 800);
    EXPECT_EQ(b.getFalseCount(), 0);

    BinaryStats c = BinaryStats::merged(a, b);

    EXPECT_NEAR(c.get(), 0.9, 1e-15);
    EXPECT_NEAR(c.getPercent(), 90.0, 1e-15);
    EXPECT_EQ(c.getTotalCount(), 1000);
    EXPECT_EQ(c.getTrueCount(), 900);
    EXPECT_EQ(c.getFalseCount(), 100);

}

TEST(RunningStats, Mean) {
    size_t const upper_limit = 2*1000*1000;
    for (size_t limit = 3; limit < upper_limit; limit = limit*4/3) {
        RunningStats s;
        for (size_t ii = 0; ii <= limit; ++ii) {
            s.push(ii);
        }
        double const dlimit = limit;
        double const squaresum = dlimit*(dlimit+1.)*(2.*dlimit+1.)/6.;
        double const sum = dlimit*(dlimit+1.)/2.;
        double const variance = ((dlimit+1.)*squaresum - sum*sum)/(dlimit*(dlimit+1.));
        EXPECT_NEAR(s.getMean(), dlimit/2, 1e-16);
        EXPECT_TRUE(RelativeNear(s.getVar(), variance, 1e-14));
    }

    for (int64_t limit = 3; limit < upper_limit; limit = limit*4/3) {
        RunningStats s;
        for (int64_t ii = -limit; ii <= limit; ++ii) {
            s.push(ii);
        }
        double const dlimit = limit;
        double const squaresum = dlimit*(dlimit+1.)*(2.*dlimit+1.)/6.;
        double const variance = squaresum/(dlimit);
        EXPECT_NEAR(s.getMean(), 0, 1e-16);
        EXPECT_TRUE(RelativeNear(s.getVar(), variance, 1e-14));
    }
}

TEST(RunningStats, Correlation) {


    double const THRESHOLD = 1e-15;
    for (long int limit = 5; limit < 1000*1000; limit *= 2) {
        {
            RunningCovariance cov;
            RunningStats stats;
            for (long int ii = -limit; ii <= limit; ++ii) {
                cov.push(ii, ii);
                stats.push(ii);
            }
            EXPECT_NEAR(cov.getMeanX(), 0, THRESHOLD);
            EXPECT_NEAR(cov.getMeanY(), 0, THRESHOLD);
            EXPECT_NEAR(cov.getCorr(), 1, THRESHOLD);

            EXPECT_NEAR(stats.getMean(), cov.getMeanX(), THRESHOLD);
            EXPECT_NEAR(stats.getVar(), cov.getVarX(), THRESHOLD);
        }
        {
            RunningCovariance cov;
            RunningStats stats;
            for (long int ii = 0; ii <= limit; ++ii) {
                cov.push(ii, -ii);
                stats.push(-ii);
            }
            double const sum = (limit * (limit + 1))/2;
            double const sum_of_squares = (sum * (2*limit + 1))/3;
            double const expected_variance = (sum_of_squares - sum*sum / (limit+1)) / limit;
            EXPECT_NEAR(cov.getMeanX(), static_cast<double>(limit)/2, THRESHOLD);
            EXPECT_NEAR(cov.getMeanY(), -static_cast<double>(limit)/2, THRESHOLD);
            EXPECT_NEAR(cov.getVarX(), expected_variance, expected_variance * THRESHOLD);
            EXPECT_NEAR(cov.getVarY(), expected_variance, expected_variance * THRESHOLD);

            EXPECT_NEAR(cov.getCorr(), -1, THRESHOLD);

            EXPECT_NEAR(stats.getMean(), cov.getMeanY(), THRESHOLD);
            EXPECT_NEAR(stats.getVar(), cov.getVarY(), THRESHOLD);
        }
        {
            RunningCovariance cov;
            RunningStats stats1, stats2;
            for (long int ii = 0; ii <= limit; ++ii) {
                cov.push(ii, limit-ii);
                stats1.push(ii);
                stats2.push(limit-ii);
            }
            double const sum = (limit * (limit + 1))/2;
            double const sum_of_squares = (sum * (2*limit + 1))/3;
            double const expected_variance = (sum_of_squares - sum*sum / (limit+1)) / limit;
            EXPECT_NEAR(cov.getMeanX(), static_cast<double>(limit)/2, THRESHOLD);
            EXPECT_NEAR(cov.getMeanY(), static_cast<double>(limit)/2, THRESHOLD);
            EXPECT_NEAR(cov.getVarX(), expected_variance, expected_variance * THRESHOLD);
            EXPECT_NEAR(cov.getVarY(), expected_variance, expected_variance * THRESHOLD);

            EXPECT_NEAR(cov.getCorr(), -1, THRESHOLD);

            EXPECT_NEAR(stats1.getMean(), static_cast<double>(limit)/2, THRESHOLD);
            EXPECT_NEAR(stats2.getMean(), static_cast<double>(limit)/2, THRESHOLD);
            EXPECT_NEAR(stats1.getVar(), expected_variance, expected_variance * THRESHOLD);
            EXPECT_NEAR(stats2.getVar(), expected_variance, expected_variance * THRESHOLD);
        }
    }

    {
        long int const limit = 1000*1000;
        RunningCovariance cov;
        RunningStats stats1, stats2;
        for (long int ii = 0; ii <= limit; ++ii) {
            cov.push(ii, limit-ii);
            stats1.push(ii);
            stats2.push(limit-ii);
        }
        double const sum = (limit * (limit + 1))/2;
        double const sum_of_squares = (sum * (2*limit + 1))/3;
        double const expected_variance = (sum_of_squares - sum*sum / (limit+1)) / limit;
        EXPECT_NEAR(cov.getMeanX(), static_cast<double>(limit)/2, static_cast<double>(limit)/2*THRESHOLD);
        EXPECT_NEAR(cov.getMeanY(), static_cast<double>(limit)/2, static_cast<double>(limit)/2*THRESHOLD);
        EXPECT_NEAR(cov.getVarX(), expected_variance, expected_variance * THRESHOLD);
        EXPECT_NEAR(cov.getVarY(), expected_variance, expected_variance * THRESHOLD);

        EXPECT_NEAR(cov.getCorr(), -1, THRESHOLD);

        EXPECT_NEAR(stats1.getMean(), static_cast<double>(limit)/2, THRESHOLD);
        EXPECT_NEAR(stats2.getMean(), static_cast<double>(limit)/2, THRESHOLD);
        EXPECT_NEAR(stats1.getVar(), expected_variance, expected_variance * THRESHOLD);
        EXPECT_NEAR(stats2.getVar(), expected_variance, expected_variance * THRESHOLD);
    }

    RunningCovariance cov;
    for (int ii = 0; ii <= 5; ++ii) {
        cov.push(ii, -ii);
    }
    cov.printInfo();
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    std::cout << "RUN_ALL_TESTS return value: " << RUN_ALL_TESTS() << std::endl;
}
