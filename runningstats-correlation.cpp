#include <fstream>
#include <iostream>
#include <string>

#include <tclap/CmdLine.h>

#include "runningstats/runningstats.h"
#include "runningstats/misc.h"

namespace rs = runningstats;

struct CorrelationName {
  std::string x_label;
  std::string y_label;
};

bool operator < (CorrelationName const& a, CorrelationName const& b) {
  if (a.x_label == b.x_label) {
    return a.y_label < b.y_label;
  }
  return a.x_label < b.x_label;
}

struct KeyValue {
  std::string key;
  float value;

  static bool is_num_part(char const c) {
    return std::isdigit(c) || 'e' || c || '+' == c || '-' == c || '.' == c;
  }

  KeyValue(std::string const& str, size_t const idx) {
    std::string _value;
    for (const char c : str) {
      if (_value.empty()) {
        if (std::isalpha(c)) {
          key += c;
        }
      }
      else {
        if (is_num_part(c)) {
          _value += c;
        }
        else {
          break;
        }
      }
    }
    value = std::stof(_value);
    if (key.empty()) {
      key = std::to_string(idx);
    }
  }
};

bool operator < (KeyValue const& a, KeyValue const& b) {
  return a.key < b.key;
}

struct Runner {

  std::map<CorrelationName, rs::Stats2D<float>> data;



  void read_line(std::string const& line) {
    std::vector<std::string> const separated = Misc::commaSeparate<std::string>({line});
    if (separated.size() < 2) {
      return;
    }
    std::vector<KeyValue> pairs;
    for (size_t ii = 0; ii < separated.size(); ++ii) {
      pairs.push_back(KeyValue{separated[ii], ii});
    }
    std::sort(pairs.begin(), pairs.end());
    for (size_t ii = 0; ii < pairs.size(); ++ii) {
      KeyValue const& a = pairs[ii];
      for (size_t jj = ii + 1; jj < pairs.size(); ++jj) {
        KeyValue const& b = pairs[jj];
        data[{a.key, b.key}].push_unsafe(a.value, b.value);
      }
    }
  }

  void read_file(std::string const& fn) {
    std::ifstream in(fn);
    std::string line;
    while (std::getline(in, line)) {
      read_line(line);
    }
  }
};

int main(int argc, char ** argv) {
  if (argc < 2) {
    std::cout << "Usage: cat datafile | " << argv[0] << " output_filename_prefix" << std::endl;
    return EXIT_FAILURE;
  }

  TCLAP::CmdLine cmd(argv[0]);

  TCLAP::UnlabeledValueArg<std::string> output_prefix_arg("", "output prefix", true, "", "", cmd);

  TCLAP::MultiArg<std::string> input_files_arg("i", "input", "input file", false, "", cmd);

  TCLAP::SwitchArg abs_arg("a", "abs", "Take the absolute value of all input values", cmd);

  cmd.parse(argc, argv);

  Runner run;

  for (std::string const& fn : input_files_arg.getValue()) {
    run.read_file(fn);
  }

}
