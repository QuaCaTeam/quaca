#include <iostream>

// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "GreensTensorFactory.h"
#include "GreensTensorPlate.h"
#include "GreensTensorVacuum.h"

// Green's tensor factory
std::shared_ptr<GreensTensor>
GreensTensorFactory::create(const std::string &input_file) {
  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // read the type of the kernel
  std::string type = root.get<std::string>("GreensTensor.type");

  // set the right pointer, show error if type is unknown
  if (type == "vacuum") {
    return std::make_shared<GreensTensorVacuum>(input_file);
  } else if (type == "plate") {
    return std::make_shared<GreensTensorPlate>(input_file);
  } else {
    std::cerr << "Error: Unknown Green's tensor type (" << type << ")!"
              << std::endl;
    exit(0);
  }
}
