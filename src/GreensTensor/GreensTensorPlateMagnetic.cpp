//
// Created by hermasim on 20/02/2020.
//
#include "GreensTensorPlateMagnetic.h"
#include <cmath>
#include <complex>

void GreensTensorPlateMagnetic::magnetic_tensor(cx_mat::fixed<3, 3> GT, Options_GreensTensor opts) {

    //Read out the relevant parameters
    double omega = opts.omega;
    double k = opts.kvec(0);
    double phi = opts.kvec(1);

    complex<double> kappa;
    // kapppa is defined to have either a purely
    // positive real part or purely negatively imaginary part
    kappa = sqrt(std::complex<double>(k_quad - omega_quad, 0.));
    kappa = std::complex<double>(std::abs(kappa.real()), -std::abs(kappa.imag()));

    complex<double> I(0,1);

    //Initialize the matrix to be diagonal
    GT.eye();
    //Set the diagonal elements
    GT *= 1. - k*cos(phi)*this->v/omega;
    GT(0,0) += k*cos(phi)*this->get_v()/omega;
    GT(1,0) += k*sin(phi)*this->get_v()/omega;
    GT(2,0) += I*kappa*this->get_v()/omega;
}
