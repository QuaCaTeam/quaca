# Friction App {docsify-ignore-all}

Besides the class [Friction](api/friction), we already implemented an application to directly calculate the noncontact friction (for an introduction see the [Tutorials](tutorials/first_calculation). If you just downloaded QuaCa, please compile it first (see [Getting Started](getting_started)). Once compiled, you find the executable `Friction` in the `bin` folder.

To calculate the noncontact friction you solely need to provide a `json` input file, as presented in [Input file API](documentation/inputfileapi). Assuming, the input file is named `todays_calculation.json` and placed in the above `data` folder, the calculation can be started by executing
```bash
quaca/bin> ./Friction --file ../data/todays_calculation.json
```
After the calculation is finished, the output will be stored in `todays_calculation.csv` at the same location as the `todays_calculation.json` file. The output contains the running variable, as for example the velocity, and the calculated friction.
