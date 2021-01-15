#include <iostream>

#include "runningstats/runningstats.h"

int main() {
    runningstats::QuantileStats<double> stats;

    double val = 0;
    while (std::cin >> val) {
        stats.push_unsafe(val);
    }

    std::cout << "# Finished reading values, 1. mean 2. stddev 3. median 4. quantile 25% 5. quantile 75% 6. mininum 7. maximum" << std::endl
              << stats.getMean() << "\t" << stats.getAccurateStddev() << "\t"
              << stats.getMedian() << "\t" << stats.getQuantile(.25) << "\t" << stats.getQuantile(.75)
              << "\t" << stats.getMin() << "\t" << stats.getMax() << std::endl;

    return EXIT_SUCCESS;
}
