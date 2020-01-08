#include "PowerSpectrumHarmOsc.h"
#include <string>
#include <armadillo>

using namespace arma;


PowerSpectrumHarmOsc::PowerSpectrumHarmOsc(std::string input_file): PowerSpectrum(input_file)
{
};

PowerSpectrumHarmOsc::PowerSpectrumHarmOsc(GreensTensor* greens_tensor, Polarizability* polarizability):PowerSpectrum(greens_tensor, polarizability)
{
};


void PowerSpectrumHarmOsc::calculate(cx_mat::fixed<3,3>& powerspectrum, double omega)
{
  cx_mat::fixed<3,3> green(fill::zeros);
  cx_mat::fixed<3,3> alpha(fill::zeros);
  Options_GreensTensor opts;
  opts.fancy_I_kv_temp = true;
  opts.omega = omega;
  this->greens_tensor->integrate_1d_k(powerspectrum, opts);
  this->polarizability->calculate( alpha, omega);
  powerspectrum = alpha*green*conj(trans(alpha));
};

