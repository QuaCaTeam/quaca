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
In this tutorial we want to compute the quantum friction force experienced by an atom moving at constant velocity $v$ and height $z_a$ above the planar surface of a semi-infinite bulk material comprosed of a simple Drude metal (as was, e.g., been done in [this publication](https://link.aps.org/doi/10.1103/PhysRevLett.123.120401)).
The computation is performed at temperature $T=\frac{1}{\beta k_\mathrm{B}}$ (see also [this publication](https://arxiv.org/abs/2110.13635) for details on the temperature-dependence of quantum friction). 
Further physical assumptions and details are explained in [this paper](http://link.aps.org/doi/10.1103/PhysRevLett.117.100402).


<div style="text-align:center">
<img class="plain" src="_media/setup.svg" class="center" width="40%">
</div>
Our setup consists of:

* a microscopic object (for example an atom), described by its polarizability $\underline{\alpha}(\omega)$
* and a vacuum-metal interface, described by the reflection coefficients $r^s(\omega,k)$ and $r^p(\omega,k)$ which are functions of the Drude permittivity and depend on the frequency $\omega$ and the wavevector $k$ of the electromagnetic field excitations.

## 2. Create the input file
Once the particle's properties - by means of the polarizability - and the material properties - by means of the permittivity, the reflection coefficients and the electric Green tensor - are defined, the friction force can be calculated. 
QuaCa exploits the inherent modularity of the friction's mathematical structure by defining a separate class for each of the input functions. 
Classes are abstract units that resemble a larger term in the formula for quantum friction such as the polarizability, the power spectrum or the Green's tensor.
Parameters for each of the keys in the class are then taken by an input file.
In this way, the program does not have to be recompiled for each calculation.

The parameters are written in  a `.json` file.
It consists of sections, which are named after the classes (e.g. `GreensTensor`, `Polarizability`) and keys which are named after the properties they represent (e.g. the velocity v, the inverse temperature beta, etc.).
Let us create an input file in the `data/` directory and call it `tutorial.json`:
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

QuaCa works with natural units, so all parameters that we enter into the input file have to be converted into natural units.

?> For more information and for a unit converter applet see the [units page](documentation/units).

Now let us dissect the input file.

First, we set the parameters for the Green's tensor, where the `type` decides on the overall geometry. In our case, this is a semi-infinite half-space given by the abbreviation `plate`.
Subsequently, we set the physical parameters of this type of Green tensor: The velocity `v`, the inverse temperature `beta` and the distance from the plate `za`. We lastly specify numerical parameters that are needed for the integration routine. You can read more about them in the [API documentation of the Green's tensor](api/greenstensor).

The next three sections specify the details of the reflection coefficients, the material's permittivity and the polarizability. In our case, we use a spatially local and homogeneous bulk material (`local bulk`), choose the Drude model (`drude`) for describing a metallic plasma with plasma frequency `omega_p` and dissipation rate `gamma`, and intend wo work with a two-level atom with resonance frequency `omega_a` and static polarizability `alpha_zero`.

In order to calculate quantum friction, we need to perform an integral over frequencies. The desired relative error of the numerical integration routine is specified in the section `Friction`.

Lastly, we tell QuaCa to calculate the frictional force for several parameters, in our case velocities. For this, we specify in the `Looper` category the type `v` together with a scale type, start/end points and number of steps.

?> For a more comprehensive explanation of all keys in a particular section, have a look at the documentation of the class, e.g. for the keys in the polarizability see [Polarizability](api/polarizability).

With our completed input file we can now start the calculation.

## 3. Run QuaCa
After [installing and building QuaCa](gettingstarted.md) you should find an executable called `Friction` in the `bin/` directory.
Change into this directory and type into the command
```bash
quaca/bin> ./Friction --file ../data/tutorial.json
```
The flag `--file` specifies the input file, where the desired path is given directly afterwards.
We choose the input file that we have created above.
Note that, since we are in the `bin/` directory, we first had to go one directory up and then to `data/tutorial.json`.
QuaCa now produces an output file in the same directory as the input file and with the same name, but of the file type `.csv`.
It contains in the first column the variable that we have looped over (which in this case is the velocity `v`) and in the second column the calculated value of the force.

## 4. Check the results
Let us now plot the data we obtained from our calculation.
We have a plot script prepared for in the `plots/` directory. 
The script reinstalls SI units and gives the frictional acceleration (force divided by mass) for a Rubidium atom as a function of velocity measured in multiples of the speed of light `c`.
Simply type
```bash
quaca/plots> python plot.py ../data/tutorial.csv
```
The result can be compared to the orange line in Fig. 2 of [this reference](https://link.aps.org/doi/10.1103/PhysRevLett.123.120401), since we have used the same parameters and worked at negligible temperatures (see input file).
