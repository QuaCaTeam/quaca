#include <iostream>
#include <sys/stat.h>

// program options
#include <boost/program_options.hpp>
namespace po = boost::program_options;

// filesystem
#include <boost/filesystem.hpp>
#include <utility>
namespace filesys = boost::filesystem;

#include "Options.h"

Options::Options(std::string parameter_file, int num_threads)
    : parameter_file(std::move(parameter_file)), num_threads(num_threads) {
  /* set according output file */
  this->output_file =
      this->parameter_file.substr(0, this->parameter_file.find_last_of('.')) +
      ".csv";

  /* check the given options */
  this->check();
}

Options::Options(int argc, char *argv[]) {
  /* Read command line options */
  try {
    // List all options and their description
    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "Help screen")(
        "file", po::value<std::string>(&(this->parameter_file)), "Input File")(
        "threads", po::value<int>(&(this->num_threads))->default_value(-1),
        "Number of parallel thrads");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // if the help option is given, show the flag description
    if (vm.count("help")) {
      std::cout << desc << std::endl;
      exit(0);
    }

  } catch (std::exception &e) {
    std::cerr << "error: " << e.what() << std::endl;
  }

  /* set according output file */
  this->output_file =
      this->parameter_file.substr(0, this->parameter_file.find_last_of('.')) +
      ".csv";

  /* check the given options */
  this->check();
}

// check the command line options
void Options::check() {
  // input file should be .json
  if (getFileExtension(parameter_file) != ".json") {
    std::cerr << "Input file must have extension .json!" << std::endl;
    exit(0);
  }

  // check if output file already exists
  if (exists(output_file)) {
    char input;
    std::cout << "Output file \"" << output_file << "\" already exists!\n"
              << "Overwrite file? [y/n]" << std::endl;
    std::cin >> input;

    if (input == 'n') {
      exit(0);
    }
  }
}

// Get File extension from File path or File Name
std::string getFileExtension(std::string filePath) {
  // Create a Path object from given string
  filesys::path pathObj(filePath);
  // Check if file name in the path object has extension
  if (pathObj.has_extension()) {
    // Fetch the extension from path object and return
    return pathObj.extension().string();
  }
  // In case of no extension return empty string
  return "";
}

// check if given file exists
bool exists(const std::string &path) {
  struct stat buffer {};
  return (stat(path.c_str(), &buffer) == 0);
}
