#include <iostream>

// ini parser
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "GreensTensorFactory.h"
#include "GreensTensorPlate.h"
#include "GreensTensorVacuum.h"

// Green's tensor factory
GreensTensor *GreensTensorFactory::create(std::string input_file) {
  // set return pointer to NULL
  GreensTensor *greenstensor = NULL;

  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // read the type of the kernel
  std::string type = root.get<std::string>("GreensTensor.type");

  // set the right pointer, show error if type is unknown
  if (type == "vacuum") {
    greenstensor = new GreensTensorVacuum(input_file);
  } else if (type == "plate") {
    greenstensor = new GreensTensorPlate(input_file);
  } else {
    std::cerr << "Error: Unknown Green's tensor type (" << type << ")!"
              << std::endl;
    exit(0);
  };

  // return memory kernel pointer
  return greenstensor;
};
