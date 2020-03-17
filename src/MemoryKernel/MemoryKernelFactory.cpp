#include "MemoryKernelFactory.h"
#include "OhmicMemoryKernel.h"
#include <iostream>

// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

// memory kernel factory
MemoryKernel *MemoryKernelFactory::create(std::string input_file,
                                          std::string section) {
  // set return pointer to NULL
  MemoryKernel *memorykernel = NULL;

  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // read the type of the kernel
  std::string type = root.get<std::string>(section + ".type");

  // set the right pointer, show error if type is unknown
  if (type == "ohmic") {
    memorykernel = new OhmicMemoryKernel(input_file, section);
  } else {
    std::cerr << "Error: Unknown Memory Kernel type (" << type << ")!"
              << std::endl;
    exit(0);
  };

  // return memory kernel pointer
  return memorykernel;
};
