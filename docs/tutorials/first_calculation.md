# First QuaCa calculation

We assume that you have QuaCa installed and are ready to go.
In this tutorial we will show you how to setup your system and how QuaCa generally works.
A calculation with QuaCa always consists of the following steps:

1. Plan the physical system
2. Create the input file
3. Run QuaCa
4. Check the results

Let us adhere to this outline and start by planing our setup.

## 1. Plan the physical system
In this tutorial we want to calculate the noncontact friction of an atom above the surface of a Drude bulk material (this 
has already been done in [this publication](https://link.aps.org/doi/10.1103/PhysRevLett.123.120401)).
Essentially we want to reproduce the results of the yellow line in FIG.2 of the aforementioned publication.

<div style="text-align:center">
<img class="plain" src="_media/setup.svg" class="center" width="40%">
</div>
Our setup consists of:

* an object (for example an atom), described by its polarizability $\underline{\alpha}(\omega)$
* a surface, described by the reflection coefficients $r^s(\omega,k)$ and $r^p(\omega,k)$
* where the object moves with constant velocity $v$ and constant distance to the surface $z_a$, and the system is evaluated at temperature $T=\frac{1}{\beta k_\mathrm{B}}$.

Further physical assumptions and details are explained in [this paper](http://link.aps.org/doi/10.1103/PhysRevLett.117.100402).
## 2. Create the input file
There are lots of things we could change in the above setup without changing the formula we have to use.
Because of this inherent modularity QuaCa reads an input file that contains all these parameters, so that the program does not have to be recompiled on each change.

The parameters are written in  a `.json` file.
It consists of sections, which are named after the classes (e.g. `GreensTensor`, `Polarizability`) and keys which are named after the properties they represent (e.g. the velocity v, the inverse temperature beta, etc.).
Classes are abstract units that resemble a larger term in the formula for quantum friction such as the polarizability, the power spectrum or the Green's tensor.
Let us now create an input file in the `data/` directory and call it `tutorial.json`:
```json
{
    "GreensTensor": {
        "type" : "plate",
        "v" : 1e-4,
        "beta" : 1e6,
        "za" : 0.025,
        "delta_cut" : 20,
        "rel_err_0" : 1e-4,
        "rel_err_1" : 1e-2
    },
    "ReflectionCoefficients": {
        "type": "local bulk"
    },
    "Permittivity": {
        "type": "drude",
        "omega_p" : 9.0,
        "gamma" : 0.035
    },
    "Polarizability": {
        "omega_a" : 1.3,
        "alpha_zero" : 6e-9
    },
    "Friction": {
        "relerr_omega" : 1e-1
    },
    "Looper": {
        "type" : "v",
        "scale" : "log",
        "start" : 1e-4,
        "end" : 1e-2,
        "steps" : 40
    }
}
```

First of all a note on units. QuaCa works with a system of measurement called natural units, so all parameters that we enter in the input file have to be converted into natural units.

?> For more information and for a unit converter applet see the [units page](documentation/units).

Now let us dissect the input file.

The first section is dedicated to the Green's tensor and contains first the type, which in our case is `plate`, the value for the physical parameters of velocity `v`, inverse temperature `beta` and the distance from the plate `za`. Furthermore we specify numerical parameters that are needed for the integration of the Green's tensor. You can read more about them in the [API documentation of the Green's tensor](api/greenstensor).

The next section is dedicated to the reflection coefficients of the surface and simply state that we want to use the reflection coefficients of a local bulk. After that we specify the permittivity of the bulk to be of type `drude`, i.e. we employ the Drude model, and specify its parameters. 

The next section defines the polarizability of the particle, which has the type `nobath` (that means the particle has no internal bath),  its parameters and the power spectrum. To calculate the quantum friction we need to perform an integral over frequencies, the desired relative error of which is specified in the section `Friction`.

Last but not least we want to tell QuaCa to calculate the frictional force for several velocities. For this we specify in the `Looper` category the type `v`together with a scale type, start/end points and number of steps.

?> For a more comprehensive explanation of all keys in a particular section, have a look the documentation of the class, e.g. for the keys in the polarizability have a look at [Polarizability](api/polarizability).

With our input file filled we can now start the calculation.

## 3. Run QuaCa
After [installing and building QuaCa](gettingstarted.md) you should find an executable called `Friction` in the `bin/` directory.
Change into this directory and type into the command line
```bash
quaca/bin> ./Friction --file ../data/tutorial.json
```
The flag `--file` specifies the input file, whose path is specified behind it.
We will chose here the input file that we have created above.
Notice that since we are in the `bin/` directory we first had to go one directory up and then to `data/tutorial.json`.
QuaCa now produces an output file in the same directory as the input file and with the same name, but of the file type `.csv`.
It contains in the first column the variable that we have looped over (which in this case is the velocity v) and in the second column the calculated value of the quantum friction.

## 4. Check the results
Let us now plot the data we obtained from our calculation and compare it to the yellowish line in [this publication](https://link.aps.org/doi/10.1103/PhysRevLett.123.120401).
We have a plot script prepared for in the `plots/` directory.
Simply type
```bash
quaca/plots> python plot.py ../data/tutorial.csv
```
and compare the plot with figure 2 of [the publication](https://link.aps.org/doi/10.1103/PhysRevLett.123.120401).
