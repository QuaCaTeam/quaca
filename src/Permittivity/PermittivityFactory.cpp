#include <iostream>

// ini parser
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "PermittivityDrude.h"
#include "PermittivityFactory.h"

// permittivity factory
Permittivity *PermittivityFactory::create(std::string input_file) {
  // set return pointer to NULL
  Permittivity *permittivity = NULL;

  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // read the type of permittivity
  std::string type = root.get<std::string>("Permittivity.type");

  // set the right pointer, show error if type is unknown
  if (type == "drude") {
    permittivity = new PermittivityDrude(input_file);
  } else {
    std::cerr << "Error: Unknown Permittivity type (" << type << ")!"
              << std::endl;
    exit(0);
  };

  // return permittivity pointer
  return permittivity;
};
