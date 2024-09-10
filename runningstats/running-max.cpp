#include "runningstats.h"

namespace runningstats {

double RunningMax::getValue() const
{
  return value;
}

void RunningMax::push_unsafe(const double val)
{
  if (std::isnan(val)) {
    return;
  }
  if (val > value) {
    value = val;
  }
}

void RunningMax::push_unsafe(const RunningStats &stats)
{
  push_unsafe(std::abs(stats.getMax()));
  push_unsafe(std::abs(stats.getMin()));
}

template<class T>
void RunningMax::push_unsafe(const QuantileStats<T> &stats)
{
  push_unsafe(std::abs(stats.getMax()));
  push_unsafe(std::abs(stats.getMin()));
}

template
void RunningMax::push_unsafe(const QuantileStats<float> &stats);

template
void RunningMax::push_unsafe(const QuantileStats<double> &stats);


} // namespace runningstats
