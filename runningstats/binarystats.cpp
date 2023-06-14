#undef NDEBUG
#include <cassert>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

#include "runningstats.h"

#include "gnuplot-iostream.h"

#include <fmt/core.h>

namespace runningstats {



void BinaryStats::push(const bool value) {
    if (value) {
        count_true++;
    }
    count_total++;
}

double BinaryStats::get() const {
    if (count_total == 0) {
        return 0;
    }
    return static_cast<double>(count_true) / count_total;
}

double BinaryStats::getPercent() const {
    return 100.0*get();
}

double BinaryStats::getFalsePercent() const {
    return 100.0*(1.0-get());
}

size_t BinaryStats::getTrueCount() const {
    return count_true;
}

size_t BinaryStats::getFalseCount() const {
    return count_total - count_true;
}

std::string BinaryStats::print() const {
    std::stringstream out;
    return fmt::format(
                "True: {:10L} ({:.4f}%), False: {:10L} ({:.4f}%), Total: {:10L}",
                count_true, getPercent(),
                getFalseCount(), getFalsePercent(),
                getTotalCount()
                );
}

void BinaryStats::clear() {
    count_total = 0;
    count_true = 0;
}

size_t BinaryStats::getTotalCount() const {
    return count_total;
}

void BinaryStats::pushTrue(size_t const num) {
    count_true += num;
    count_total += num;
}

void BinaryStats::pushFalse(size_t const num) {
    count_total += num;
}

BinaryStats BinaryStats::merged(BinaryStats const& a, BinaryStats const& b) {
    BinaryStats result(a);
    result.pushTrue(b.getTrueCount());
    result.pushFalse(b.getFalseCount());

    return result;
}

} // namespace runningstats
