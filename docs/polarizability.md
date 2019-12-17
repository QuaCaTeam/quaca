# Polarizability

With Spike it is possible to simulate several neuron types.
They're all specified in the input file (format `.ini`) with the appropriate parameters.

QuaCa offers several polarizabilities.
All are of the form 
$ \underline{\alpha} = \alpha_0 \left( \omega_a^2 - \omega - i \omega \mu(\omega) \right)^{-1}$

## PolarizabilityBath
Implements a polarizability with an internal bath.
Obeys the equation
$ \underline{\alpha} = \alpha_0 \left( \omega_a^2 - \omega - i \omega \mu(\omega) \right)^{-1}$.


## PolarizabilityNoBath
Implements a polarizability without an internal bath.
Obeys the equation
$ \underline{\alpha} = \alpha_0 \left( \omega_a^2 - \omega \right)^{-1}$.
