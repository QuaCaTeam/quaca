#ifndef POWERSPECTRUMHARMOSC
#define POWERSPECTRUMHARMOSC

#include "PowerSpectrum.h"
#include <string>

class PowerSpectrumHarmOsc : public PowerSpectrum {
public:
  PowerSpectrumHarmOsc(GreensTensor *greens_tensor,
                       Polarizability *polarizability);
  PowerSpectrumHarmOsc(std::string input_file);

  // calculate the power spectrum for a fixed value of the frequency
  void calculate(cx_mat::fixed<3, 3> &powerspectrum, double omega);
};

#endif
