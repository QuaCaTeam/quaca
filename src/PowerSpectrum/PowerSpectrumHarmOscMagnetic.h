//
// Created by hermasim on 03/03/2020.
//

#ifndef QUACA_POWERSPECTRUMHARMOSCMAGNETIC_H
#define QUACA_POWERSPECTRUMHARMOSCMAGNETIC_H

#include "PowerSpectrumHarmOsc.h"
#include "PowerSpectrum.h"
#include "../GreensTensor/GreensTensorPlateMagnetic.h"


class PowerSpectrumHarmOscMagnetic: public PowerSpectrumHarmOsc {
public:
  GreensTensorPlateMagnetic* greens_magnetic;
    //Constructor with initalization list
    PowerSpectrumHarmOscMagnetic(GreensTensor *greens_tensor,
    Polarizability *polarizability);
    //Constructor with ini-file
    PowerSpectrumHarmOscMagnetic(std::string input_file);

    // calculate the power spectrum for a fixed value of the frequency
    void calculate(cx_mat::fixed<3, 3> &powerspectrum,
                   Options_PowerSpectrum opts);
};


#endif //QUACA_POWERSPECTRUMHARMOSCMAGNETIC_H
