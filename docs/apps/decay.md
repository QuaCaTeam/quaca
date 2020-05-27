# Decay Rate App {docsify-ignore-all}

To calculate the decay rate give via
$$
\Gamma(\omega) = \alpha_0\omega_a^2\mathrm{Tr}\Big[\frac{1}{\underline{\alpha}(\omega)}\cdot\underline{\alpha}_\Im(\omega)\cdot\frac{1}{\underline{\alpha}^\dagger(\omega)}\Big] 
$$
you can use the `Decay rate` app, which is already implemented in QuaCa. If you just downloaded QuaCa, please compile it first (see [Getting Started](getting_started)). Once compiled, you find the executable `Decay` in the `bin` folder.

To calculate the decay rate you solely need to provide a `json` input file, as presented in [Input file API](documentation/inputfileapi). Assuming, the input file is named `todays_calculation.json` and placed in the above `data` folder, the calculation can be started by executing
```bash
quaca/bin> ./Decay --file ../data/todays_calculation.json
```
After the calculation is finished, the output will be stored in `todays_calculation.csv` at the same location as the `todays_calculation.json` file. The output contains the running variable $\omega$, and the calculated friction.
