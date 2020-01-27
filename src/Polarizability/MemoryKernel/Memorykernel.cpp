// ini parser
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "MemoryKernel.h"

MemoryKernel::MemoryKernel(){};

MemoryKernel::MemoryKernel(std::string input_file) {
  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // read parameters
  this->type = root.get<std::string>("MemoryKernel.type");
}
