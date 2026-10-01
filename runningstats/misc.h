#ifndef MISC_H
#define MISC_H

#include <string>
#include <vector>

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

    static void trim(std::string &s);

    template<class T>
    static std::vector<T> commaSeparate(std::vector<std::string> const &args)
    {
      std::vector<T> result;
      for (std::string const &s : args) {
        std::string current;
        for (char c : s) {
          if (c == ',' || c == ';' || c == ' ') {
            trim(current);
            if (!current.empty()) {
              std::stringstream stream(current);
              T val;
              stream >> val;
              result.push_back(val);
              current = "";
            }
          }
          else {
            current += c;
          }
        }
        trim(current);
        if (!current.empty()) {
          std::stringstream stream(current);
          T val;
          stream >> val;
          result.push_back(val);
        }
      }
      return result;
    }
};

#endif // MISC_H
