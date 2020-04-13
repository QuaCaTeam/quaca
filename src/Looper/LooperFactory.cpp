#include <iostream>

// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "LooperFactory.h"
#include "LooperV.h"
#include "LooperZa.h"

std::shared_ptr<Looper> LooperFactory::create(const std::string& input_file) {
  // Create a root
  // Load the json file in this ptree
  pt::ptree root;
  pt::read_json(input_file, root);

  // read the type of the kernel
  std::string type = root.get<std::string>("Looper.type");

  // set the right pointer, show error if type is unknown
  if (type == "v") {
    return std::make_shared<LooperV>(input_file);
  } else if (type == "za") {
    return std::make_shared<LooperZa>(input_file);
  } else {
    std::cerr << "Error: Unknown Looper type (" << type << ")!" << std::endl;
    exit(0);
  }
}
