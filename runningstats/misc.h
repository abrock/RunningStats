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
};

#endif // MISC_H
