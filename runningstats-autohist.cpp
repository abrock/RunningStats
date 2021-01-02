#include <iostream>

#include "runningstats/runningstats.h"

template<class T>
double getBinSize(runningstats::QuantileStats<T> & stats) {
    double const iqr = stats.getQuantile(.75) - stats.getQuantile(.25);
    return 2*iqr/(std::cbrt(double(stats.getCount())));
}

int main(int argc, char ** argv) {
    if (argc < 2) {
        std::cout << "Usage: cat datafile | " << argv[0] << " output_filename_prefix" << std::endl;
        return EXIT_FAILURE;
    }

    runningstats::QuantileStats<double> stats;

    double val = 0;
    while (std::cin >> val) {
        stats.push_unsafe(val);
    }

    std::cout << "# Finished reading values, 1. mean 2. stddev 3. median 4. quantile 25% 5. quantile 75%" << std::endl
              << stats.getMean() << "\t" << stats.getAccurateStddev() << "\t"
              << stats.getMedian() << "\t" << stats.getQuantile(.25) << "\t" << stats.getQuantile(.75) << std::endl;
    stats.plotHistAndCDF(argv[1], getBinSize(stats));

    return EXIT_SUCCESS;
}
