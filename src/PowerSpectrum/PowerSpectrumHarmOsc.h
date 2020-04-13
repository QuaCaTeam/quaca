#ifndef POWERSPECTRUMHARMOSC_H
#define POWERSPECTRUMHARMOSC_H

#include "PowerSpectrum.h"
#include <string>

class PowerSpectrumHarmOsc : public PowerSpectrum {
public:
  // Constructor with initalization list
  PowerSpectrumHarmOsc(std::shared_ptr<GreensTensor> greens_tensor,
                       std::shared_ptr<Polarizability> polarizability);
  // Constructor with json-file
  PowerSpectrumHarmOsc(const std::string& input_file);

  // calculate the power spectrum for a fixed value of the frequency
  void calculate(cx_mat::fixed<3, 3> &powerspectrum,
                 Options_PowerSpectrum opts);
  bool has_bath;
};

#endif // POWERSPECTRUMHARMOSC_H
