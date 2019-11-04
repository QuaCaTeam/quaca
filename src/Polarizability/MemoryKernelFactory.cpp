#include <iostream>

// ini parser
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
namespace pt = boost::property_tree;

#include "MemoryKernelFactory.h"
#include "OhmicMemoryKernel.h"

MemoryKernel * MemoryKernelFactory::create(std::string input_file)
{
    // set return pointer to NULL
    MemoryKernel *memorykernel = NULL;

    // Create a root
    pt::ptree root;

    // Load the ini file in this ptree
    pt::read_ini(input_file, root);

    // read simulation data into simulation variables
    std::string type = root.get<std::string>("MemoryKernel.type");

    // set the right pointer, show error if type is unknown
    if ( type == "ohmic" )
    {
      memorykernel = new OhmicMemoryKernel(input_file);
    }
    else
    {
        std::cerr << "Error: Unknown Memory Kernel type (" << type << ")!" <<std::endl;
        exit(0);
    };

    // return memory kernel pointer
    return memorykernel;
};
