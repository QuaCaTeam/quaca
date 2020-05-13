# First QuaCa calculation {docsify-ignore-all}

We assume that you have QuaCa installed and are ready to go.
In this tutorial we will show you how to setup your system and how QuaCa generally works.
A calculation always consists of the following steps:

1. Plan the physical system
2. Populate the input file
3. Run QuaCa
4. Check the results

Let us adhere to this outline and start by planing our setup.

## 1. Plan the physical system
In this tutorial we want to calculate the non-contact friction of an atom above the surface of a Drude bulk material, which has already been done in [this publication](https://link.aps.org/doi/10.1103/PhysRevLett.123.120401).
So essentially we want to reproduce the results of the yellowish line in FIG.2 of the aforementioned publication.

<div style="text-align:center">
<img class="plain" src="_media/setup.svg" class="center" width="40%">
</div>
Our setup consists of:

* an object (for example an atom), described by its polarizability $\underline{\alpha}(\omega)$
* a surface, described by the reflection coefficients $r^s(\omega,k)$ and $r^p(\omega,k)$
* where the object moves with constant velocity $v$ and constant distance to the surface $z_a$, and the system is evaluated at temperature $T=\frac{1}{\beta k_\mathrm{B}}$.

Further physical assumptions and subtilities are explained in [this paper](http://link.aps.org/doi/10.1103/PhysRevLett.117.100402).
## 2. Populate the input file
There are lots of things we could change in the above setup without changing the formula we have to use.
Because of this inherent modularity QuaCa reads an input file that contains all these parameters, so that the program does not have to be recompiled on each change.

The parameters are written in  a `.json` file.
It consists of sections, which are named after the classes (e.g. GreensTensor, Polarizability) and keys which are named after the properties they represent (e.g. the velocity v, the inverse temperature beta, etc.).
Classes are abstract units that resemble a larger term in the formula for quantum friction such as the polarizability, the power spectrum or the Green's tensor.
For our situation the section for the Green's tensor would look like this:
```json
{
    "GreensTensor":
    {
        "type": "plate",
        "v" : 5e-4,
        "beta": 1e-4
    }
}
```

?> To see which keys are available in a certain section, look into the [API documentation](api) of the relevant class.

Let us now create an input file in the `data/` directory and call it `tutorial.json`.
Let us now list all classes that we have to define in the input file.
We need a Green's tensor for which we in turn need reflection coefficients and a permittivity.
```json
{
    "GreensTensor":{
        ...
    },
    "ReflectionCoefficients":{
        ...
    },
    "Permittivity":{
        ...
    }
}
```

We also need a polarizability and a power spectrum, so our input file looks like this
```json
{
    "GreensTensor":{
        ...
    },
    "ReflectionCoefficients":{
        ...
    },
    "Permittivity":{
        ...
    },
    "Polarizability":{

    },
    "PowerSpectrum":{

    }
}
```
Let us now specify each class according to our description above and add the keys (i.e. the variable names) to the appropriate section
```json
{
    "GreensTensor":{
        "type" : "plate",
        "v" : ,
        "beta" : ,
        "za" : ,
        "delta_cut" : ,
        "rel_err_0" : ,
        "rel_err_1" : 
    },
    "ReflectionCoefficients":{
        "type": "local bulk"
    },
    "Permittivity":{
        "type": "drude",
        "omega_p" : ,
        "gamma" : 
    },
    "Polarizability":{
        "type" : "nobath",
        "omega_a" : ,
        "alpha_zero" :
    },
    "PowerSpectrum":{
        "type" : "harmonic oscillator"
    }
}
```

You might notice some parameters that have nothing to do with the physical system we described above.
The parameters `delta_cut`, `rel_err_0` and `rel_err_1` are numerical parameters that specify the cutoff for the integrations up to infinity and the relative errors of the integration respectively.
For more information have a look at the documentation of the [GreensTensor class](api/greenstensor).

We are now ready to enter some values into our input file, but wait.
What are the units that QuaCa requires?
QuaCa works with a system of measurement called natural units.
For more information and for a unit converter applet see the [units page](documentation/units).
If we convert the units for our system appropriately the input file `tutorial.json` should now look like this
```json
{
    "GreensTensor":{
        "type" : "plate",
        "v" : 1e-4,
        "beta" : 1e6,
        "za" : 0.01,
        "delta_cut" : 20,
        "rel_err_0" : 1e-4,
        "rel_err_1" : 1e-2
    },
    "ReflectionCoefficients":{
        "type": "local bulk"
    },
    "Permittivity":{
        "type": "drude",
        "omega_p" : 9.0,
        "gamma" : 0.1
    },
    "Polarizability":{
        "type" : "nobath",
        "omega_a" : 1.3,
        "alpha_zero" : 6e-9
    },
    "PowerSpectrum":{
        "type" : "harmonic oscillator"
    }
}
```

Now that we have entered all values into our input file we are almost ready to start our first calculation.
The only thing left to specify is that we want to calculate the quantum friction over a range of velocities.
This can be specified by defining a *Looper* at the end of the input file.
Let us define the Looper at the end of our input file, so our final version looks like this
```json
{
    "GreensTensor":{
        "type" : "plate",
        "v" : 1e-4,
        "beta" : 1e6,
        "za" : 0.01,
        "delta_cut" : 20,
        "rel_err_0" : 1e-4,
        "rel_err_1" : 1e-2
    },
    "ReflectionCoefficients":{
        "type": "local bulk"
    },
    "Permittivity":{
        "type": "drude",
        "omega_p" : 9.0,
        "gamma" : 0.1
    },
    "Polarizability":{
        "type" : "nobath",
        "omega_a" : 1.3,
        "alpha_zero" : 6e-9
    },
    "PowerSpectrum":{
        "type" : "harmonic oscillator"
    },
    "Friction":{
        "relerr_omega" : 1e-1
    },
    "Looper":{
        "type" : "v",
        "scale" : "log",
        "start" : 1e-4,
        "end" : 1e-2,
        "steps" : 40
    }
}
```
With our input file filled we can now finally start the calculation.

## 3. Run QuaCa
After installing and building QuaCa you should find an executable called `Friction` in the `bin/` directory.
Change into this directory and type into the command line
```bash
quaca/bin> ./Friction --file ../data/tutorial.json
```
The flag `--file` specifies the input file, whose path is specified behind it.
We will chose here the input file that we have populated above.
Notice that since we are in the `bin/` directory we first had to go one directory up and then to `data/tutorial.json`.
QuaCa now produces an output file in the same directory as the input file and with the same name, but of the file type `.csv`.
It contains in the first column the variable that we have looped over (which in this case is the velocity v) and in the second column the calculated value of the quantum friction.

The calculation might take a while, with the above paramters and depending on your computer around 10-25sec per value (on the PCs of the developers this calculation ran in 10-14min). If you do not want to wait that long, reduce the number of steps in the input file.

## 4. Check the results
Let us now plot the data we obtained from our calculation and compare it to the yellowish line in [this publication](https://link.aps.org/doi/10.1103/PhysRevLett.123.120401).
We have a plot script prepared for in the `plots/` directory.
Simply type
```bash
quaca/plots> python plot.py ../data/tutorial.csv
```
and compare the plot with figure 2 of [the publication](https://link.aps.org/doi/10.1103/PhysRevLett.123.120401).
