#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

double getMedian(std::string const& filename) {
    std::ifstream in(filename);
    std::string line;
    double last_val = 0;
    while (std::getline(in, line)) {
        std::stringstream _line(line);
        double value = 0;
        double cdf = 0;
        if (!_line) {
            continue;
        }
        _line >> value;
        if (!_line) {
            continue;
        }
        _line >> cdf;
        last_val = value;
        if (cdf >= 0.5) {
            return value;
        }
    }
    return last_val;
}

double pow10(double const x) {
    double res = std::pow(10.0, x);
    res = std::round(res*1000)/1000;
    return res;
}

int main(int argc, char ** argv) {

    for (size_t ii = 1; ii < size_t(argc); ++ii) {
        std::cout << argv[ii] << " " << pow10(getMedian(argv[ii])) << std::endl;
    }

    return 0;
}
