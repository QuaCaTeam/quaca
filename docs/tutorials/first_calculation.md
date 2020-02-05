!> TODO: fill in values and finish writing

# First QuaCa calculation

We assume that you have QuaCa installed and are ready to go.
In this tutorial we will show you how to setup your system and how QuaCa generally works.
A calculation always consists of the following steps:

1. Plan the physical system
2. Populate the input file
3. Run QuaCa
4. Check the results

Let us adhere to this outline and start by planing our setup.

## 1. Plan the physical system
In this tutorial we want to calculate *describe calculation* , which has already been done in *give reference to publication*.
So essentially we want to reproduce the plot *blabla on page blublu*.

Our setup consists of:

* a particle, described by
* a material filling the lower halfspace with

Furthermore we assume that ...

## 2. Populate the input file
There are lots of things we could change in the above setup without changing the formula we have to use.
For example, imagine we want to calculate the exact same situation as above, but with a *blabla* atom.
Because of this inherent modularity QuaCa reads an input file that contains all these parameters, so that the program does not have to be recompiled on each change.

The parameters are written in  a `.ini` file.
It consists of sections, which are named after the classes (e.g. GreensTensor, Polarizability) and keys which are named after the properties they represent (e.g. the velocity v, the inverse temperature beta, etc.).
Classes are abstract units that resemble a larger term in the formula for quantum friction such as the polarizability, the power spectrum or the Green's tensor.
For our situation the section for the Green's tensor would look like this:
```ini
[GreensTensor]
type = plate
v = 5e-4
beta = 1e-4
```

?> To see which keys are available in a certain section, look into the [API documentation](api) of the relevant class.

Let us now create an input file in the `data/` directory and call it `tutorial.ini`.
Let us now list all classes that we have to define in the input file
```ini
[Permittivity]

[GreensTensor]

[Polarizability]
```

According to our setup these classes need the following parameters
```ini
[Permittivity]
type = drude
omega_p =
gamma =

[GreensTensor]
type = plate
v =
beta =
z_a =

[Polarizability]
type = nobath
omega_a =
alpha_zero =
```

We are now ready to enter some values into our input file, but wait.
What are the units that QuaCa requires?
QuaCa works with a system of measurement called natural units.
For more information and for a unit converter applet see the [units page](documentation/units).
If we convert the units for our system appropriately the input file `tutorial.ini` should now look like this
```ini
[Permittivity]
type = drude
omega_p = 9.0
gamma = 35e-3

[GreensTensor]
type = plate
v = 5e-4
beta = 1.16e13
z_a = 5e-9

[Polarizability]
type = bath
omega_a = 1.3
alpha_zero = 6e-9
```

We have created our input file and are now ready to start our first calculation.

## 3. Run QuaCa
After installing and building QuaCa you should find an executable called `QuaCa` in the `bin/` directory.
Change into this directory and type into the command line
```bash
quaca/bin> ./QuaCa --file ../data/tutorial.ini
```
The flag `--file` specifies the input file, whose path is specified behind it.
Notice that since we are in the `bin/` directory we first had to go one directory up and then to `data/tutorial.ini`.
QuaCa now produces a file in the same directory as the input file and with the same name, but of the file type `.csv`.
It contains in the first column that variable that we have looped over (which in this case is XXX) and in the second column the calculated value of the quantum friction.

## 4. Check the results
Let us now plot the data we obtained from our calculation and compare it to *plot in paper*.
