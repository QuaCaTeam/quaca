//
// Created by hermasim on 03/03/2020.
//

#include "PowerSpectrumHarmOscMagnetic.h"
#include "../Polarizability/PolarizabilityNoBathMagnetic.h"
#include <armadillo>
#include <string>

using namespace arma;

//Constructor with ini-file
PowerSpectrumHarmOscMagnetic::PowerSpectrumHarmOscMagnetic(std::string input_file)
        : PowerSpectrumHarmOsc(input_file) {
	greens_magnetic = dynamic_cast<GreensTensorPlateMagnetic* >(greens_tensor);
	};

//Constructor with initialization list
PowerSpectrumHarmOscMagnetic::PowerSpectrumHarmOscMagnetic(GreensTensor *greens_tensor,
                                           Polarizability *polarizability)
        : PowerSpectrumHarmOsc(greens_tensor, polarizability){
	greens_magnetic = dynamic_cast<GreensTensorPlateMagnetic* >(greens_tensor);
	};

//Compute the power spectrum for a given frequency \omega
void PowerSpectrumHarmOscMagnetic::calculate(cx_mat::fixed<3, 3> &powerspectrum,
                                     Options_PowerSpectrum opts) {
    //Read out variables double omega = opts.omega;
    double omega = opts.omega;

    //Compute the full spectrum
    if (opts.spectrum == FULL) {
        //Initialize tensor storing the Green's tensor and setting the integration
        //options for the Green's tensor
        cx_mat::fixed<3, 3> green(fill::zeros);
        Options_GreensTensorMagnetic opts_g;
        opts_g.fancy_complex = IM;
        opts_g.BE = IM;
        opts_g.EB = IM;
        opts_g.BB = IM;
        opts_g.weight_function = TEMP;
        opts_g.omega = omega;
        opts_g.class_pt = this->greens_magnetic;

        // Compute the Green's tensor
        this->greens_magnetic->integrate_k(green, opts_g);

        //Initialize tensor storing the polarizability and seting the integration
        //options for the polarizability
        cx_mat::fixed<3, 3> alpha(fill::zeros);
        Options_Polarizability opts_alpha;
        opts_alpha.omega = omega;

        // Compute the polarizability
        this->polarizability->calculate_tensor(alpha, opts_alpha);

        // Combine the Green's tensor and the polarizability, see eq. [3.9] in
        // Marty's PhD thesis
        powerspectrum = 1. / M_PI * alpha * green * trans(alpha);

        //check wether polarizability has an internal bath
        if(has_bath) {
            std::cerr << "This feature is not implemented at this moment. Please check the documentation.";
            exit(0);
        }
    }

    //Compute only the non-LTE contributions to the power spectrum
    if (opts.spectrum == NON_LTE_ONLY) {
        std::cerr << "This feature is not implemented for the moment. Please check the documentation.";
        exit(0);
    }
};
