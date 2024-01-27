#include "runningstats.h"

namespace runningstats {

MeanFromHist::MeanFromHist() {}

MeanFromHist &MeanFromHist::operator =(const MeanFromHist &other) {
    sum = other.sum;
    n = other.n;
    return *this;
}

void MeanFromHist::push(const double value, const size_t count) {
    if (!std::isfinite(value)) {
        return;
    }
    std::lock_guard guard(push_mutex);
    push_unsafe(value, count);
}

void MeanFromHist::push_unsafe(const double value, const size_t count) {
    if (!std::isfinite(value)) {
        return;
    }
    sum += value*count;
    n += count;
    n_bins++;
}

double MeanFromHist::getMean() const {
    if (0 == n) {
        return 0;
    }
    return sum/n;
}

size_t MeanFromHist::getCount() const {
    return n;
}

size_t MeanFromHist::getBinCount() const {
    return n_bins;
}

MeanFromHist::MeanFromHist(const MeanFromHist &rhs) {
    *this = rhs;
}

} // namespace runningstats
