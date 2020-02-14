#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <iostream>
#include <omp.h>
namespace pt = boost::property_tree;

#include "ProgressBar.hpp"
#include "Quaca.h"

int main(int argc, char *argv[]) {
  // get file from command line
  Options opts(argc, argv);
  std::string parameters = opts.get_parameter_file();

  // define looper
  Looper *looper = LooperFactory::create(parameters);

  // define progressbar
  ProgressBar progbar(looper->get_steps_total(), 70);
  progbar.display();
  //Create a paralle region given threads given by the --threads flag
  //we have to create the parallel region already here to ensure,
  //that any thread creates their own instance of quantum_friction
  #pragma omp parallel num_threads(opts.get_num_threads())
  {
  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(parameters, root);

  double relerr_omega = root.get<double>("Friction.relerr_omega");

  // define needed quantities
  Polarizability *polarizability = PolarizabilityFactory::create(parameters);
  PowerSpectrumHarmOsc *powerspectrum = new PowerSpectrumHarmOsc(
      polarizability->get_greens_tensor(), polarizability);
  Friction *quant_friction =
      new Friction(polarizability->get_greens_tensor(), polarizability,
                          powerspectrum, relerr_omega);


  //array to store all the computed values of the loop
  double friction_data[looper->get_steps_total()];

  //Parallelize the for-loop of the given looper
    #pragma omp for 
    for (int i = 0; i < looper->get_steps_total(); i++) {
      std::cout << "Computing step with thread number: " << omp_get_thread_num();
      std::cout << " of " << omp_get_num_threads() << std::endl;
      friction_data[i] = looper->calculate_value(i, quantum_friction);

      std::cout << "Step: " << i << " Friction force: " << friction_data[i] << endl;
      ++progbar;
      progbar.display();
  }
  }

  //Store calculated data
  
  // define output file
  std::ofstream file;
  file.open(opts.get_output_file());

  //write the data to the file
  double step;
  for (int i = 0; i < looper->get_steps_total(); i++) {
    step = looper->get_step(i);
    file << step << "," << friction_data[i] << "\n";
  }

  // close file
  file.close();

  // close progress bar
  progbar.done();

  return 0;
};
