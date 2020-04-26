#ifndef POWERSPECTRUMHARMOSC_H
#define POWERSPECTRUMHARMOSC_H

#include "PowerSpectrum.h"
#include <string>

class PowerSpectrumHarmOsc : public PowerSpectrum {
public:
  bool has_bath;

  // Constructor with initalization list
  PowerSpectrumHarmOsc(std::shared_ptr<GreensTensor> greens_tensor,
                       std::shared_ptr<Polarizability> polarizability);
  // Constructor with json-file
  PowerSpectrumHarmOsc(const std::string &input_file);

  // calculate the power spectrum for a fixed value of the frequency
  void calculate(double omega, cx_mat::fixed<3, 3> &powerspectrum,
                 Spectrum_Options spectrum) const;
};

#endif // POWERSPECTRUMHARMOSC_H
