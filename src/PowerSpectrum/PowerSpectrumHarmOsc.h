#ifndef POWERSPECTRUMHARMOSC_H
#define POWERSPECTRUMHARMOSC_H

#include "PowerSpectrum.h"
#include <string>

class PowerSpectrumHarmOsc : public PowerSpectrum {
public:
  PowerSpectrumHarmOsc(GreensTensor *greens_tensor,
                       Polarizability *polarizability);
  PowerSpectrumHarmOsc(std::string input_file);

  // calculate the power spectrum for a fixed value of the frequency
  void calculate(cx_mat::fixed<3, 3> &powerspectrum,
                 Options_PowerSpectrum opts);
};

#endif // POWERSPECTRUMHARMOSC_H
