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

TEST(Image2D, plot) {
    Image2D<double> test(1, 0.025);
    for (double xx = -10; xx <= 10; ++xx) {
        for (double yy = -10; yy <= 10; yy += 0.025) {
            test[xx][yy] = yy;
        }
    }
    test.plot("x,y->y", HistConfig().setXLabel("x").setYLabel("y").setTitle("f(x,y) := y"));
}

TEST(Image2D, contours) {
    double const step = 0.0125;
    double const limit = 3;
    Image2D<double> test(step, step);
    for (double xx = -limit; xx <= limit; xx += step) {
        for (double yy = -limit; yy <= limit; yy += step) {
            test[xx][yy] = std::exp(-(xx*xx+yy*yy));
        }
    }
    test.plot("gaussian", HistConfig().setXLabel("x").setYLabel("y").setTitle("f(x,y) := e^(-(x^2+y^2))").addContour(.5, ".5", "ffffff"));
}

TEST(Image2D, lines) {
    double const step = 0.0125;
    double const limit = 3;
    Image2D<double> test(step, step);
    for (double xx = -limit; xx <= limit; xx += step) {
        for (double yy = -limit; yy <= limit; yy += step) {
            test[xx][yy] = std::exp(-(xx*xx+yy*yy));
        }
    }
    LineSegment l1, l2;
    for (double angle = 0; angle <= 1; angle += .02) {
        l1.addPt(2*std::cos(angle), 2*std::sin(angle));
        l2.addPt(2*std::cos(M_PI + angle), 2*std::sin(M_PI + angle));
    }
    Line l;
    l.color = "white";
    l.addSegment(l1);
    l.addSegment(l2);
    HistConfig conf;
    conf
            .setXLabel("x")
            .setYLabel("y")
            .setTitle("f(x,y) := e^(-(x^2+y^2))")
            .addLine(l)
            .setColorMap(ColorMaps::viridis_clipped());
    HistConfig copy = conf.clone();
    test.plot("gaussian+lines", copy);
}

TEST(Image2D, plot_float) {
    Image2D<float> test(1, 1);
    for (double xx = 0; xx <= 100; ++xx) {
        for (double yy = 0; yy <= 100; yy++) {
            test[xx][yy] = yy;
        }
    }
    test.plot("float x,y->y", HistConfig().setXLabel("x").setYLabel("y").setTitle("f(x,y) := y"));
}

TEST(Image2D, extractor) {
    Image2D<QuantileStats<float> > test(1, 0.025);
    for (double xx = -10; xx <= 10; ++xx) {
        for (double yy = -10; yy <= 10; yy += 0.025) {
            test[xx][yy].push_unsafe(yy);
        }
    }
    test.plot("x,y->y", HistConfig().setXLabel("x").setYLabel("y").setTitle("f(x,y) := y").extractMedian());
}

TEST(Image2d, merged) {
    QuantileStats<float> stats;
    Image2D<QuantileStats<float> > img(.1,.1);
    std::mt19937_64 engine(0xCAFECAFE);
    std::uniform_real_distribution<float> dist(-1,1);
    for (size_t ii = 0; ii < 10'000; ++ii) {
        float const val = dist(engine);
        float const x = dist(engine);
        float const y = dist(engine);
        img[x][y].push_unsafe(val);
        stats.push_unsafe(val);
    }
    QuantileStats<float> cmp_stats = img.merged<float>();
    double const threshold = 1e-7;
    ASSERT_NEAR(stats.getMedian(), cmp_stats.getMedian(), threshold);
    ASSERT_NEAR(stats.getMean(), cmp_stats.getMean(), threshold);
    ASSERT_NEAR(stats.getStddev(), cmp_stats.getStddev(), threshold);
    ASSERT_NEAR(stats.getVar(), cmp_stats.getVar(), threshold);
    for (size_t ii = 0; ii <= 200; ++ii) {
        double const quantile = double(ii)/200;
        ASSERT_NEAR(stats.getQuantile(quantile), cmp_stats.getQuantile(quantile), threshold) << "Quantile: " << quantile;
    }

}

TEST(Image2D, push_unsafe) {
    Image2D<std::vector<QuantileStats<float> > > test(1,1);
    test.push_unsafe(0,0,{0,1,2,3,4});
    test.push_unsafe(1,0,{0,1,2,3,4});
    test.push_unsafe(0,1,{0,1,2,3,4});
    test.push_unsafe(1,1,{0,1,2,3,4});
    std::cout << "Mean values: " << std::endl;
    test.data2file(std::cout, HistConfig().extractMean());
    std::cout << "Stddev: " << std::endl;
    test.data2file(std::cout, HistConfig().extractStddev());
}

TEST(Image2D, minmax) {
    {
        Image2D<double> test(1, 1);
        test[0][0] = 1;
        EXPECT_EQ(test.min_x, 0);
        EXPECT_EQ(test.max_x, 0);
        EXPECT_EQ(test.min_y, 0);
        EXPECT_EQ(test.max_y, 0);
    }
    {
        Image2D<double> test(1, 1);
        for (double ii = -10; ii <= 10; ++ii) {
            test[ii][0] = 1;
            EXPECT_EQ(test.min_x, -10);
            EXPECT_EQ(test.max_x, ii);
            EXPECT_EQ(test.min_y, 0);
            EXPECT_EQ(test.max_y, 0);
        }
        for (double ii = -10; ii <= 10; ++ii) {
            test[0][ii] = 1;
            EXPECT_EQ(test.min_x, -10);
            EXPECT_EQ(test.max_x, 10);
            EXPECT_EQ(test.min_y, -10);
            EXPECT_EQ(test.max_y, std::max(0.0, ii));
        }
    }
}

TEST(Image2D, colormaps) {
    Image2D<double> test(1, 1);
    // Fill with a simple gradient
    for (int xx = -50; xx < 250; ++xx) {
        double value = double(xx)/200;
        value = std::min<double>(1, std::max<double>(0, value));
        for (int yy = -10; yy < 20; ++yy) {
            test[xx][yy] = value;
        }
    }
    HistConfig conf;
    conf.setTitle("Simple gradient experiment");
    for (std::string const& map :
    {
         ColorMaps::blue_red(),
         ColorMaps::blue_red_clipped(),
         ColorMaps::blue_red_2(),
         ColorMaps::viridis(),
         ColorMaps::viridis_clipped(),
    }
         )
    {
    conf.setColorMap(map);
        test.plot("gradient-" + map, conf.clone().setTitle("Gradient with colormap " + map));
        test.plot("gradient-unclipped-" + map, conf.clone().setTitle("Unclipped gradient with colormap " + map).setMinMaxCB(-.01, 1.01));
    }
}

TEST(Image1D, plot) {
    {
        Image1D<RunningStats> test(1);
        for (int ii = 0; ii <= 20; ++ii) {
            test[ii].push_unsafe(ii);
            test[ii].push_unsafe(ii-ii);
            test[ii].push_unsafe(ii+ii);
        }
        test.plot("Image1D-with-errorbars", HistConfig().extractMeanAndStddev());
    }
    {
        Image1D<QuantileStats<float> > test(1);
        for (double ii = 0; ii <= 20; ++ii) {
            for (size_t jj = 0; jj < 10; ++jj) {
                test[ii].push_unsafe(ii);
            }
            for (size_t jj = 0; jj < 10; ++jj) {
                test[ii].push_unsafe(ii-ii/4);
            }
            for (size_t jj = 0; jj < 10; ++jj) {
                test[ii].push_unsafe(ii+ii/2);
            }
        }
        test.plot("Image1D-with-quantiles", HistConfig().extractMedianAndIQR());
    }
}


TEST(Image1D, minmax) {
    {
        Image1D<double> test(1);
        test[0] = 1;
        EXPECT_EQ(test.min_val, 0);
        EXPECT_EQ(test.max_val, 0);
        test[15] = 1;
        EXPECT_EQ(test.min_val, 0);
        EXPECT_EQ(test.max_val, 15);
        test[-2.2] = 1;
        EXPECT_EQ(test.min_val, -2);
        EXPECT_EQ(test.max_val, 15);
    }
    {
        Image1D<double> test(1);
        for (double ii = -10; ii <= 10; ++ii) {
            test[ii] = 1;
            EXPECT_EQ(test.min_val, -10);
            EXPECT_EQ(test.max_val, ii);
        }
    }
    {
        Image1D<double> test(1);
        for (double ii = 10; ii >= 10; --ii) {
            test[ii] = 1;
            EXPECT_EQ(test.max_val, 10);
            EXPECT_EQ(test.min_val, ii);
        }
    }
    {
        Image1D<double> test(5.5);
        for (double ii = 10; ii >= 10; --ii) {
            test[ii] = 1;
            EXPECT_EQ(test.max_val, std::round(10.0/5.5)*5.5);
            EXPECT_EQ(test.min_val, std::round(ii/5.5)*5.5);
        }
    }
}
