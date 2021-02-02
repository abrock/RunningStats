#include <iostream>
#include <random>

#include <gtest/gtest.h>

#include "runningstats/runningstats.h"

using namespace runningstats;

static std::random_device dev;
static std::mt19937 rng(dev());


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

TEST(BinaryStats, empty) {
    BinaryStats bin;
    bin.get();
    EXPECT_EQ(0, bin.getFalseCount());
    EXPECT_EQ(0, bin.getTotalCount());
    EXPECT_EQ(0, bin.getTrueCount());
}

TEST(RunningStats, empty) {
    RunningStats s;
    EXPECT_EQ(0, s.getCount());
    EXPECT_EQ(0, s.getMean());
    EXPECT_EQ(0, s.getStddev());
    EXPECT_EQ(0, s.getVar());
}

TEST(QuantileStats, empty) {
    QuantileStats<float> s;
    EXPECT_EQ(0, s.getCount());
    EXPECT_EQ(0, s.getMean());
    EXPECT_EQ(0, s.getStddev());
    EXPECT_EQ(0, s.getVar());
    EXPECT_EQ(0, s.getAccurateStddev());
    EXPECT_EQ(0, s.getAccurateVariance());
    EXPECT_EQ(std::vector<float>(), s.getData());
    s.getHistogram(1);
    s.FreedmanDiaconisBinSize();
    s.plotHistAndCDF("tmp", s.FreedmanDiaconisBinSize());
    EXPECT_EQ(0, s.getInverseQuantile(0));
    EXPECT_EQ(0, s.getLogMean());
    EXPECT_EQ(0, s.getLogStddev());
    EXPECT_EQ(0, s.getLogVar());
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

    for (int64_t limit = 3; limit < int64_t(upper_limit); limit = limit*4/3) {
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

TEST(RunningCorrelation, minmax) {
    runningstats::RunningCovariance c;

    EXPECT_EQ(0, c.getMinX());
    EXPECT_EQ(0, c.getMinY());
    EXPECT_EQ(0, c.getMaxX());
    EXPECT_EQ(0, c.getMaxY());

    c.push(0,0);

    EXPECT_EQ(0, c.getMinX());
    EXPECT_EQ(0, c.getMinY());
    EXPECT_EQ(0, c.getMaxX());
    EXPECT_EQ(0, c.getMaxY());

    c.push(0,1);

    EXPECT_EQ(0, c.getMinX());
    EXPECT_EQ(0, c.getMinY());
    EXPECT_EQ(0, c.getMaxX());
    EXPECT_EQ(1, c.getMaxY());

    c.push(1,1);

    EXPECT_EQ(0, c.getMinX());
    EXPECT_EQ(0, c.getMinY());
    EXPECT_EQ(1, c.getMaxX());
    EXPECT_EQ(1, c.getMaxY());

    c.push(1,1);

    EXPECT_EQ(0, c.getMinX());
    EXPECT_EQ(0, c.getMinY());
    EXPECT_EQ(1, c.getMaxX());
    EXPECT_EQ(1, c.getMaxY());

    c.push(-1,1);

    EXPECT_EQ(-1, c.getMinX());
    EXPECT_EQ(0, c.getMinY());
    EXPECT_EQ(1, c.getMaxX());
    EXPECT_EQ(1, c.getMaxY());

    c.push(-1,-1);

    EXPECT_EQ(-1, c.getMinX());
    EXPECT_EQ(-1, c.getMinY());
    EXPECT_EQ(1, c.getMaxX());
    EXPECT_EQ(1, c.getMaxY());

    c.push(-2,-2);

    EXPECT_EQ(-2, c.getMinX());
    EXPECT_EQ(-2, c.getMinY());
    EXPECT_EQ(1, c.getMaxX());
    EXPECT_EQ(1, c.getMaxY());

    c.push(-3,3);

    EXPECT_EQ(-3, c.getMinX());
    EXPECT_EQ(-2, c.getMinY());
    EXPECT_EQ(1, c.getMaxX());
    EXPECT_EQ(3, c.getMaxY());

    c.push(3,4);

    EXPECT_EQ(-3, c.getMinX());
    EXPECT_EQ(-2, c.getMinY());
    EXPECT_EQ(3, c.getMaxX());
    EXPECT_EQ(4, c.getMaxY());
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
        int64_t const limit = 1000*1000;
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

TEST(Histogram, all) {
    std::vector<double> const bin_sizes = {.25, .3, .5, 1., 1.4, 2.8, 5.6};
    for (const double bin_size : bin_sizes) {
        Histogram h(bin_size);
        RunningStats compare;
        for (size_t ii = 0; ii < 10; ++ii) {
            h.push(double(ii)*bin_size);
            compare.push(double(ii)*bin_size);
        }
        EXPECT_NEAR(h.getMean(), compare.getMean(), 1e-16);
        EXPECT_NEAR(h.getVar(), compare.getVar(), 1e-16);
        EXPECT_NEAR(h.getMin(), compare.getMin(), 1e-16);
        EXPECT_NEAR(h.getMax(), compare.getMax(), 1e-16);
        EXPECT_NEAR(h.getCount(), compare.getCount(), 1e-16);
        EXPECT_NEAR(h.getStddev(), compare.getStddev(), 1e-16);
        EXPECT_NEAR(h.getLogMean(), compare.getLogMean(), 1e-16);
        EXPECT_NEAR(h.getLogStddev(), compare.getLogStddev(), 1e-16);
        {
            auto hist = h.getAbsoluteHist();
            size_t counter = 0;
            for (auto & entry : hist) {
                EXPECT_NEAR(double(counter)*bin_size, entry.first, 1e-16);
                EXPECT_NEAR(1, entry.second, 1e-16);

                ++counter;
            }
        }
        {
            auto hist = h.getRelativeHist();
            size_t counter = 0;
            for (auto & entry : hist) {
                EXPECT_NEAR(double(counter)*bin_size, entry.first, 1e-16);
                EXPECT_NEAR(1.0 / (bin_size*10), entry.second, 1e-16);

                ++counter;
            }
        }
    }
}

TEST(plot, plotHist) {
    runningstats::Histogram h(0.1);
    std::normal_distribution<double> dist;
    for (size_t ii = 0; ii < 500; ++ii) {
        h.push(dist(rng));
    }
    h.plotHist("test-plot-hist", false);
}

TEST(plot, QuantileStats_plotHist) {
    runningstats::QuantileStats<double> h;
    std::normal_distribution<double> dist;
    for (size_t ii = 0; ii < 5'000; ++ii) {
        h.push(dist(rng));
    }
    h.plotHist("test-plot-quantile-hist", 0.1, false);
    h.plotCDF("test-plot-quantile-cdf");
    h.plotCDF("test-plot-quantile-cdf-absolute", HistConfig().setAbsolute());
    HistConfig conf;
    conf.setXLabel("random value").setYLabel("density").setTitle("Testing HistConf");
    h.plotHist("test-plot-quantile-hist-conf", h.FreedmanDiaconisBinSize(), conf);

    for (size_t ii = 0; ii < 95*1000; ++ii) {
        h.push(dist(rng));
    }
    h.plotHist("test-plot-quantile-hist-100k", h.FreedmanDiaconisBinSize(), false);
#pragma omp parallel sections
    {
#pragma omp section
        h.plotCDF("test-plot-quantile-cdf-100k");
#pragma omp section
        h.plotCDF("test-plot-quantile-cdf-100k-absolute", HistConfig().setAbsolute());
    }
    conf.setXLabel("random value").setYLabel("density").setTitle("Testing HistConf");
    h.plotHist("test-plot-quantile-hist-conf-100k", h.FreedmanDiaconisBinSize(), conf);

    for (size_t ii = 0; ii < 900*1000; ++ii) {
        h.push(dist(rng));
    }
    h.plotHist("test-plot-quantile-hist-1M", h.FreedmanDiaconisBinSize(), false);
    h.plotCDF("test-plot-quantile-cdf-1M");
    conf.setXLabel("random value").setYLabel("density").setTitle("Testing HistConf");
    h.plotHist("test-plot-quantile-hist-conf-1M", h.FreedmanDiaconisBinSize(), conf);
}

TEST(plot, QuantileStats_plotHist_range) {
    runningstats::QuantileStats<double> h;
    std::normal_distribution<double> dist;
    for (size_t ii = 0; ii < 5000; ++ii) {
        h.push(dist(rng));
    }
    h.plotHist("test-plot-quantile-minmax-hist", 0.1, false);
    h.plotCDF("test-plot-quantile-minmax-cdf");
    HistConfig conf;
    conf.setIgnoreAmount(.1);
    conf.setXLabel("random value").setYLabel("density").setTitle("Testing HistConf");
    h.plotHist("test-plot-quantile-minmax-hist-conf", h.FreedmanDiaconisBinSize(), conf);

    for (size_t ii = 0; ii < 95*1000; ++ii) {
        h.push(dist(rng));
    }
    h.plotHist("test-plot-quantile-minmax-hist-100k", h.FreedmanDiaconisBinSize(), false);
    h.plotCDF("test-plot-quantile-minmax-cdf-100k", conf);
    conf.setXLabel("random value").setYLabel("density").setTitle("Testing HistConf");
    h.plotHist("test-plot-quantile-minmax-hist-conf-100k", h.FreedmanDiaconisBinSize(), conf);

    for (size_t ii = 0; ii < 900*1000; ++ii) {
        h.push(dist(rng));
    }
    h.plotHist("test-plot-quantile-minmax-hist-1M", h.FreedmanDiaconisBinSize(), false);
    h.plotCDF("test-plot-quantile-minmax-cdf-1M", conf);
    conf.setXLabel("random value").setYLabel("density").setTitle("Testing HistConf");
    h.plotHist("test-plot-quantile-minmax-hist-conf-1M", h.FreedmanDiaconisBinSize(), conf);
}

TEST(plot, QuantileStats_plotHistAndCDF) {
    runningstats::QuantileStats<double> h;
    std::normal_distribution<double> dist;

    for (size_t ii = 0; ii < 1000*1000; ++ii) {
        h.push(dist(rng));
    }
    h.plotHistAndCDF("test-plot-hist-and-cdf-1M", h.FreedmanDiaconisBinSize(), HistConfig().setDataLabel("N(0,1)"));
}

TEST(std_vector, capacity) {
    for (size_t ii = 0; ii < 20; ++ii) {
        std::vector<size_t> vec;
        vec.reserve(ii);
        std::cout << "Reserved: " << ii << ", capacity: " << vec.capacity() << std::endl;
    }
    for (size_t ii = 0; ii < 20; ++ii) {
        std::vector<size_t> vec;
        for (size_t jj = 0; jj < ii; ++jj) {
            vec.push_back(jj);
        }
        std::cout << "Pushed: " << ii << ", capacity: " << vec.capacity() << std::endl;
    }
    for (size_t ii = 0; ii < 20; ++ii) {
        std::vector<size_t> vec(ii, 0);
        std::vector<int> copy(vec.begin(), vec.end());
        std::cout << "Initialized: " << ii << ", capacity: " << vec.capacity() << ", copy: " << copy.capacity() << std::endl;
    }
}

TEST(Ellipses, Stats2D_eillpises) {
    runningstats::Ellipses ellipses;
    std::normal_distribution<double> dist;

    for (size_t ii = 0; ii < 20; ++ii) {
        runningstats::Stats2D<double> stats;
        for (size_t jj = 0; jj < 1000; ++jj) {
            stats.push_unsafe(dist(rng), dist(rng));
        }
        stats.getQuantileEllipse(ellipses, 0.5);
    }
    ellipses.plot("ellipses-stats2d", HistConfig()
                  .setTitle("50% ellipses of 1k normal distributions")
                  .setXLabel("X label")
                  .setYLabel("Y label")
                  .setFixedRatio());
}

TEST(Ellipses, plot) {
    runningstats::Ellipses ellipses;
    for (size_t xx = 1; xx < 10; ++xx) {
        ellipses.push({xx, 0, 2*xx, 10});
    }
    ellipses.plot("ellipses", HistConfig()
                  .setTitle("Ellipses")
                  .setXLabel("X label")
                  .setYLabel("Y label"));
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    std::cout << "RUN_ALL_TESTS return value: " << RUN_ALL_TESTS() << std::endl;
}
