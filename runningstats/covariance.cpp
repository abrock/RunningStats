#undef DNDEBUG
#include <cassert>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

#include "runningstats.h"

#include "gnuplot-iostream.h"

namespace runningstats {

void RunningCovariance::push_unsafe(double x, double y)
{
    n++;
    if (n == 1) {
        meanX = maxX = minX = x;
        meanY = maxY = minY = y;
    }
    else {
        double const newMeanX = meanX + (x - meanX) / n;
        double const newMeanY = meanY + (y - meanY) / n;

        varSumX = varSumX + (x - meanX) * (x - newMeanX);
        varSumY = varSumY + (y - meanY) * (y - newMeanY);

        double const newCovarSum = covarSum + (x - meanX) * (y - meanY) * (n-1) / n;

        meanX = newMeanX;
        meanY = newMeanY;
        covarSum = newCovarSum;

        minX = std::min(minX, x);
        minY = std::min(minY, y);
        maxX = std::max(maxX, x);
        maxY = std::max(maxY, y);
    }
}

RunningCovariance::RunningCovariance() {}

RunningCovariance::RunningCovariance(const RunningCovariance &rhs) {
    *this = rhs;
}

RunningCovariance &RunningCovariance::operator=(const RunningCovariance &other) {
    n = other.n;
    minX = other.minX;
    maxX = other.maxX;
    minY = other.minY;
    maxY = other.maxY;
    varSumX = other.varSumX;
    varSumY = other.varSumY;
    covarSum = other.covarSum;
    return *this;
}

void RunningCovariance::push(double x, double y) {
    std::lock_guard<std::mutex> const guard(push_mutex);
    push_unsafe(x, y);
}

double RunningCovariance::getMeanX() {
    return meanX;
}

double RunningCovariance::getVarX() const {
    return n <= 1 ? 0 : varSumX/(n - 1);
}

double RunningCovariance::getMeanY() {
    return meanY;
}

double RunningCovariance::getVarY() const {
    return n <= 1 ? 0 : varSumY/(n - 1);
}

double RunningCovariance::getStddevX() const {
    return std::sqrt(getVarX());
}

double RunningCovariance::getStddevY() const {
    return std::sqrt(getVarY());
}

double RunningCovariance::getCoVar() {
    return covarSum / (n-1);
}

double RunningCovariance::getCorr() {
    return getCoVar() / std::sqrt(getVarX() * getVarY());
}

double RunningCovariance::getMinX() const {
    return minX;
}

double RunningCovariance::getMaxX() const {
    return maxX;
}

double RunningCovariance::getMinY() const {
    return minY;
}

double RunningCovariance::getMaxY() const {
    return maxY;
}

void RunningCovariance::printInfo() {
    std::cout.precision(15);
    std::cout << std::scientific
              << "Variables: " << std::endl
              << "meanX: " << meanX << std::endl
              << "meanY: " << meanY << std::endl
              << "varSumX: " << varSumX << std::endl
              << "varSumY: " << varSumY << std::endl
              << "covarSum: " << covarSum << std::endl
              << "getVarX(): " << getVarX() << std::endl
              << "getVarY(): " << getVarY() << std::endl
              << "getCoVar(): " << getCoVar() << std::endl
              << "getCorr(): " << getCorr() << std::endl
              << std::endl;
}

} // namespace runningstats
