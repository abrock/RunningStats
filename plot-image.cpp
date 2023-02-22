#include <iostream>

#include <opencv2/highgui.hpp>

#include <tclap/CmdLine.h>

#include "runningstats/runningstats.h"
namespace rs = runningstats;

int const grid_cells = 400;

int main(int argc, char ** argv) {
    TCLAP::CmdLine cmd("Plot pixel values of greyscale images");

    TCLAP::UnlabeledMultiArg<std::string> input_arg("", "Input images", true, "Image filenames", cmd);

    cmd.parse(argc, argv);

#pragma omp parallel for ordered schedule(dynamic, 1)
    for (size_t ii = 0; ii < input_arg.getValue().size(); ++ii) {
        std::string const& input = input_arg.getValue()[ii];
        cv::Mat_<float> const img = cv::imread(input, cv::IMREAD_ANYDEPTH | cv::IMREAD_GRAYSCALE);

        rs::Image2D<rs::QuantileStats<float> > out(double(img.cols)/grid_cells, double(img.rows)/grid_cells);
        rs::RunningStats stats;

        for (int yy = 0; yy < img.rows; ++yy) {
            for (int xx = 0; xx < img.cols; ++xx) {
                out[xx][yy].push_unsafe(img(yy,xx));
                stats.push_unsafe(img(yy,xx));
            }
        }

        rs::HistConfig conf;
        conf
                .extractMean()
                .setFixedRatio()
                .setFlipY()
                .setMinMaxCB(stats.getMin(), stats.getMax());

        out.plot(input, conf);
#pragma omp ordered
        {
            std::cout << "Stats for " << input << ": " << std::endl
                      << stats.printBoth() << std::endl;
        }
    }


    return EXIT_SUCCESS;
}
