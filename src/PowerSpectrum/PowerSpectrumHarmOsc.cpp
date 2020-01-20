#include "PowerSpectrumHarmOsc.h"
#include <armadillo>
#include <string>

using namespace arma;

PowerSpectrumHarmOsc::PowerSpectrumHarmOsc(std::string input_file)
    : PowerSpectrum(input_file){};

PowerSpectrumHarmOsc::PowerSpectrumHarmOsc(GreensTensor *greens_tensor,
                                           Polarizability *polarizability)
    : PowerSpectrum(greens_tensor, polarizability){};

void PowerSpectrumHarmOsc::calculate(cx_mat::fixed<3, 3> &powerspectrum,
                                     double omega) {
  cx_mat::fixed<3, 3> green(fill::zeros);
  Options_GreensTensor opts_g;
  opts_g.fancy_I_kv_temp = true;
  opts_g.omega = omega;
  opts_g.class_pt = this->greens_tensor;
  this->greens_tensor->integrate_1d_k(green, opts_g);

  cx_mat::fixed<3, 3> alpha(fill::zeros);
  Options_Polarizability opts_alpha;
  opts_alpha.omega = omega;
  this->polarizability->calculate_tensor(alpha, opts_alpha);

  powerspectrum = alpha * green * conj(trans(alpha));
};
