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

    std::cout << "Finished reading values" << std::endl;
    stats.plotHistAndCDF(argv[1], getBinSize(stats));

    return EXIT_SUCCESS;
}
