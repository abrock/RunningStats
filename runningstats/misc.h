#ifndef MISC_H
#define MISC_H

#include <string>

struct Misc
{
    /**
     * @brief make_target_dir makes sure the directory for a given filename exists.
     * @param filename
     */
    static void make_target_dir(std::string const& filename);

    /**
     * @brief range2string creates a range-string for usage in xrange, yrange, cbrange etc.
     * If one of the values is not finite it will be represented by an empty string, letting gnuplot automatically decide the value.
     * @param min
     * @param max
     * @return
     */
    static std::string range2string(double const min, double const max);

    static std::string replace_empty(std::string const& val, std::string const& replacement);
};

#endif // MISC_H
