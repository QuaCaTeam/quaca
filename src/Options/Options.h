#ifndef OPTIONS_H
#define OPTIONS_H

#include <string>

// Class for command line options
class Options {
private:
  std::string parameter_file; // path to parameter file
  std::string output_file;    // path to output file

public:
  // constructor
  Options(std::string parameter_file);
  Options(int argc, char *argv[]);

  // check input options
  void check();

  // getter function
  std::string get_parameter_file() { return this->parameter_file; };
  std::string get_output_file() { return this->output_file; };
};

// get extension of a file
std::string getFileExtension(std::string filePath);

// check if file exists
bool exists(const std::string &path);

#endif // OPTIONS_H
