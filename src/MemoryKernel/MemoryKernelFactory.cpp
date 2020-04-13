#include "MemoryKernelFactory.h"
#include "OhmicMemoryKernel.h"
#include <iostream>

// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

// memory kernel factory
std::shared_ptr<MemoryKernel>
MemoryKernelFactory::create(const std::string &input_file,
                            const std::string &section) {
  // Create a root
  // Load the json file in this ptree
  pt::ptree root;
  pt::read_json(input_file, root);

  // read the type of the kernel
  std::string type = root.get<std::string>(section + ".type");

  // set the right pointer, show error if type is unknown
  if (type == "ohmic") {
    return std::make_shared<OhmicMemoryKernel>(input_file, section);
  } else {
    std::cerr << "Error: Unknown Memory Kernel type (" << type << ")!"
              << std::endl;
    exit(0);
  }
}
