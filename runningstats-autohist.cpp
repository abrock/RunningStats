#include <iostream>
#include <fstream>

#include "runningstats/runningstats.h"
namespace rs = runningstats;

#include <tclap/CmdLine.h>

double conditional_abs(double const value, bool const do_abs) {
    return do_abs ? std::abs(value) : value;
}

int main(int argc, char ** argv) {
    if (argc < 2) {
        std::cout << "Usage: cat datafile | " << argv[0] << " output_filename_prefix" << std::endl;
        return EXIT_FAILURE;
    }

    TCLAP::CmdLine cmd(argv[0]);

    TCLAP::UnlabeledValueArg<std::string> output_prefix_arg("", "output prefix", true, "", "", cmd);

    TCLAP::MultiArg<std::string> input_files_arg("i", "input", "input file", false, "", cmd);

    TCLAP::SwitchArg abs_arg("a", "abs", "Take the absolute value of all input values", cmd);

    TCLAP::SwitchArg log_x_arg("", "log-x", "Use logscale for the x-axis", cmd);

    TCLAP::SwitchArg log_y_arg("", "log-y", "Use logscale for the y-axis", cmd);

    cmd.parse(argc, argv);

    bool const do_abs = abs_arg.getValue();

    runningstats::QuantileStats<double> stats;



    for (std::string const& file : input_files_arg.getValue()) {
        double val = 0;
        if ("--" == file) {
            while (std::cin >> val) {
                stats.push_unsafe(conditional_abs(val, do_abs));
            }
        }
        else {
            std::ifstream in(file);
            while (in >> val) {
                stats.push_unsafe(conditional_abs(val, do_abs));
            }
        }
    }

    std::string const prefix = output_prefix_arg.getValue();

    rs::HistConfig conf;
    if (log_x_arg.getValue()) {
        conf.setLogX();
    }
    if (log_y_arg.getValue()) {
        conf.setLogY();
    }

    std::cout << "# Finished reading values, 1. mean 2. stddev 3. median 4. quantile 25% 5. quantile 75%" << std::endl
              << stats.getMean() << "\t" << stats.getAccurateStddev() << "\t"
              << stats.getMedian() << "\t" << stats.getQuantile(.25) << "\t" << stats.getQuantile(.75) << std::endl;
    stats.plotHistAndCDF(prefix, conf);
    stats.plotHistAndCDF(prefix + "-0.01", conf.clone().setIgnoreAmount(.01));
    stats.plotHistAndCDF(prefix + "-0.05", conf.clone().setIgnoreAmount(.05));
    stats.plotHistAndCDF(prefix + "-0.10", conf.clone().setIgnoreAmount(.10));
    stats.plotHistAndCDF(prefix + "-0.20", conf.clone().setIgnoreAmount(.20));

    return EXIT_SUCCESS;
}
