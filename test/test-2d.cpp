#include <iostream>
#include <random>

#include <gtest/gtest.h>

#include "runningstats/runningstats.h"

using namespace runningstats;

#include "randutils.hpp"

static std::random_device dev;
static std::mt19937 engine(dev());

TEST(Plot2D, Histogram2D) {
    std::normal_distribution<double> dist;

    runningstats::Histogram2D hist(0.1, 0.2);

    for (size_t ii = 0; ii < 1000*1000; ++ii) {
        hist.push_unsafe(.2*dist(engine), dist(engine));
    }

    hist.plotHist("normal-2d", false);

}

TEST(Plot2D, Stats2D) {
    std::normal_distribution<double> dist;

    runningstats::Stats2D<float> stats;

    for (size_t ii = 0; ii < 1000*1000; ++ii) {
        stats.push_unsafe(.1*dist(engine), dist(engine));
    }

    runningstats::Histogram2D hist = stats.getHistogram2D(stats.FreedmanDiaconisBinSize());
    hist.plotHist("stats-normal-2d", false);

}

TEST(Plot2D, LaplaceHist) {
    std::normal_distribution<double> dist;

    runningstats::Stats2D<float> stats;

    randutils::mt19937_rng rng;

    for (size_t ii = 0; ii < 1'000'000; ++ii) {
        float const v1 = rng.variate<float, std::uniform_real_distribution>(-1,1);
        stats.push_unsafe(
                    v1,
                    v1+rng.variate<float, std::uniform_real_distribution>(-.1,.1)
                    );
    }
    std::pair<double, double> bin = {.01, .01};
    stats.plotHist("plot2d-laplace", bin,
                   HistConfig()
                   .addExtractor(HistConfig::Extract::Median, 0, "black", "Median")
                   .addExtractor(HistConfig::Extract::Quantile, .75, "white", "IQR")
                   .addExtractor(HistConfig::Extract::Quantile, .25, "white", "")
                   );
    stats.getHistogram2Dfixed(bin).plotHistPm3D("plot2d-laplace-pm3d");
    stats.getHistogram2Dfixed(bin)
            .plotHistPm3D("plot2d-laplace-pm3d-normalized",
                          HistConfig()
                          .setNormalizeX()
                          .addExtractors({{HistConfig::Extract::Mean, 0},
                                         {HistConfig::Extract::Quantile, 0.75},
                                         {HistConfig::Extract::Quantile, 0.25}}));

}

TEST(Plot2D, Stats2Dnormalized) {
    std::normal_distribution<double> dist;

    runningstats::Stats2D<float> stats;

    for (size_t ii = 0; ii < 1000*1000; ++ii) {
        stats.push_unsafe(.1*dist(engine), dist(engine));
    }
    stats.plotHist("stats-normal-2d-normalized", stats.FreedmanDiaconisBinSize(),
                   HistConfig().setNormalizeX().setLogCB());
}

TEST(Plot2D, Stats2Dlimited) {
    std::normal_distribution<double> dist;

    runningstats::Stats2D<float> stats;

    for (size_t ii = 0; ii < 1000*1000; ++ii) {
        stats.push_unsafe(.1*dist(engine), dist(engine));
    }

    stats.plotHist("stats-normal-2d-limited", stats.FreedmanDiaconisBinSize(), HistConfig().setMaxBins(10,10));
}

TEST(Plot2D, Stats2Dfixed) {
    std::normal_distribution<double> dist;

    runningstats::Stats2D<float> stats;

    for (size_t ii = 0; ii < 1000*1000; ++ii) {
        stats.push_unsafe(.1*dist(engine), dist(engine));
    }

    std::pair<double, double> bin = stats.FreedmanDiaconisBinSize();
    HistConfig conf;
    conf.setTitle("2D normal distribution, F-D bin width").setXLabel("N(0,0.1²)").setYLabel("N(0,1²)").setLogCB();
    stats.plotHist("stats-fixed-normal-2d", bin, conf);

    conf.setIgnoreAmount(.001);
    conf.setTitle("2D normal distribution, F-D bin width").setXLabel("N(0,0.1²)").setYLabel("N(0,1²)").setLogCB();
    stats.plotHist("stats-fixed-normal-2d-ignore-0.001", bin, conf);
}

TEST(Plot2D, Stats2Dfixed_large) {
    std::normal_distribution<double> dist;

    runningstats::Stats2D<float> stats;

    for (size_t ii = 0; ii < 10*1000*1000; ++ii) {
        stats.push_unsafe(.1*dist(engine), dist(engine));
    }

    std::pair<double, double> bin = stats.FreedmanDiaconisBinSize();
    HistConfig conf;
    conf.setTitle("2D normal distribution, F-D bin width").setXLabel("N(0,0.1²)").setYLabel("N(0,1²)").setLogCB();
    stats.plotHist("stats-fixed-normal-2d-large", bin, conf);

    conf.setIgnoreAmount(.001);
    conf.setTitle("2D normal distribution, F-D bin width").setXLabel("N(0,0.1²)").setYLabel("N(0,1²)").setLogCB();
    stats.plotHist("stats-fixed-normal-2d-large-ignore-0.001", bin, conf);
}

TEST(Plot2D, Stats2Dfixed_degenerate) {
    std::normal_distribution<double> dist;

    runningstats::Stats2D<float> stats;

    for (size_t ii = 0; ii < 10*1000; ++ii) {
        stats.push_unsafe(0, dist(engine));
    }

    std::pair<double, double> bin = stats.FreedmanDiaconisBinSize();
    HistConfig conf;
    conf.setTitle("2D normal distribution, F-D bin width").setXLabel("0").setYLabel("N(0,1²)").setLogCB();
    stats.getHistogram2Dfixed(bin).plotHist("stats-fixed-degenerate-2d", conf);
    conf.setTitle("2D normal distribution, F-D/2 bin width");
    stats.getHistogram2Dfixed({bin.first/2, bin.second/2}).plotHist("stats-fixed-degenerate-2d-half", conf);

}

TEST(Stats2D, merge) {
    Stats2D<float> a, b, c;
    for (size_t ii = 0; ii < 10; ++ii) {
        a.push_unsafe(0,1);
        b.push_unsafe(1,0);
    }
    c.push(a);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    std::cout << "RUN_ALL_TESTS return value: " << RUN_ALL_TESTS() << std::endl;
}
