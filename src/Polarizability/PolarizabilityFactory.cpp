#include <iostream>

// ini parser
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
namespace pt = boost::property_tree;

#include "PolarizabilityFactory.h"
#include "PolarizabilityBath.h"
#include "PolarizabilityNoBath.h"

// Green's tensor factory
Polarizability * PolarizabilityFactory::create(std::string input_file)
{
    // set return pointer to NULL
    Polarizability *polarizability = NULL;

    // Create a root
    pt::ptree root;

    // Load the ini file in this ptree
    pt::read_ini(input_file, root);

    // read the type of the kernel
    std::string type = root.get<std::string>("Polarizability.type");

    // set the right pointer, show error if type is unknown
    if ( type == "bath" )
    {
      polarizability = new PolarizabilityBath(input_file);
    }
    else if ( type == "nobath" )
    {
      polarizability = new PolarizabilityNoBath(input_file);
    }
    else
    {
        std::cerr << "Error: Unknown Polarizability type (" << type << ")!" <<std::endl;
        exit(0);
    };

    // return memory kernel pointer
    return polarizability;
};
