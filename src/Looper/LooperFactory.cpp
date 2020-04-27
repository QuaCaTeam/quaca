#include <iostream>

// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "LooperFactory.h"
#include "LooperV.h"
#include "LooperZa.h"

Looper *LooperFactory::create(std::string input_file) {
  // set return pointer to NULL
  Looper *looper = NULL;

  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // read the type of the kernel
  std::string type = root.get<std::string>("Looper.type");

  // set the right pointer, show error if type is unknown
  if (type == "v") {
    looper = new LooperV(input_file);
  } else if (type == "za") {
    looper = new LooperZa(input_file);
  } else if (type == "omega") {
    looper = new LooperOmega(input_file);
  } else {
    std::cerr << "Error: Unknown Looper type (" << type << ")!" << std::endl;
    exit(0);
  };

  // return looper pointer
  return looper;
};
