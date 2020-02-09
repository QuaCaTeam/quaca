!> TODO: Fill in values and finish writing

# Write your own main
In the first tutorial you calculated the non-contact friction of an atom above the surface of a Drude bulk material using an input file and the executable `Friction`.
QuaCa can, however, also be used like a library and you can write your own executables!
In this tutorial we will redo the last tutorial, but this time we will implement everything in our own main file.
You will learn how to:

1. Register a new executable.
2. Write the main function.

## 1. Register a new executable
Let us look at the subfolders of `quaca/`.
```
    quaca/
    ├── app
    ├── bin
    ├── build
    ├── data
    ├── docs
    ├── include
    ├── lib
    ├── notes
    ├── src
    └── test
```

You will see amongst them a folder called `app/` and a folder called `src/`.
The folder `src/` contains all the code for the Green's tensors, polarizabilities and all other quantities that are needed for the calculation of quantum friction.
This folder does, however, *not* contain a main function, which is needed to build an executable (i.e. a program we can actually execute).
These main functions are inside the `app/` directory.
Look inside the `app/` folder and you will find a subfolder called `QuaCa` and a file called `CMakeLists.txt`.

Each folder inside the `app/` folder should contain a `.cpp` file containing a main file, so each folder corresponds to an executable program.
In the last tutorial you have already worked with the executable `QuaCa`, so let us now create a new one.
First create a new folder inside the `app/` folder and call it `Tutorial`.
Next go to the file `CMakeLists.txt` inside the `app/` directory and add the following line to the end of the file
```cmake
add_subdirectory("Tutorial")
```
This tells CMake to consider this folder when building the executables.
Next create a file inside the `Tutorial/` folder called `tutorial.cpp` and write the following inside it.
```cpp
#include "Quaca.h"
#include <iostream>

int main(int argc, char *argv[])
{
  std::cout << "The tutorial executable works." << std::endl;
  return 0;
};
```
Now we need to tell CMake how to build this executable.
Create a file called `CMakeLists.txt` inside the `Tutorial/` folder and write the following in it
```cmake
# add executable
add_executable(Tutorial
  ${CMAKE_CURRENT_SOURCE_DIR}/tutorial.cpp
  )
target_include_directories(Tutorial PRIVATE ../include)

# link libraries
target_link_libraries(Tutorial
  quaca
  )
```
This tells CMake to add an executable called `Tutorial`, that the file containing the main function is called `tutorial.cpp` and that this executable should be linked to the library `quaca`.

With this setup ready we should test whether our tutorial works.
Go to the `quaca/build/` directory and type into the command line
```bash
quaca/build> cmake ..
```
You have to do this again, because we created new files and CMake needs to know about them!
The build the project typing
```bash
quaca/build> make
```
You should now see your executable in the `quaca/bin/` directory, so let us execute it and we see
```bash
quaca/bin> ./Tutorial
The tutorial executable works!
quaca/bin>
```
Good job!
Let us now proceed by filling our main with more meaningfull code.

## 2. Write the main function
First let us define the permittivity.
We can do so by typing
```cpp
#include "Quaca.h"
#include <iostream>

int main(int argc, char *argv[])
{
  // define permittivity
  double omega_p = 9.0;
  double gamma = 0.1;
  PermittivityDrude permittivity(omega_p, gamma);

  return 0;
};
```
You have now created an object of the type PermittivityDrude called `permittivity` with the appropriate parameters.
Let us proceed by defining our Green's tensor and put this permittivity inside it.
Type below the permittivity definition
```cpp
// define the Green's tensor
double v = 1e-4;
double beta = 1e6;
double z_a = 0.01;

double delta_cut = 20;
vec::fixed<2> rel_err = {1E-4, 1E-2};

GreensTensorPlate greens_tensor(v, z_a, beta, &permittivity);
```
We have now defined the Green's tensor and have it the permittivity which we defined above it.

This process has to be repeated with the polarizability
```cpp
// define polarizability
double omega_a = 1.3;
double alpha_zero = 6e-9;
PolarizabilityBath polarizability(omega_a, alpha_zero, &greens_tensor);
```
With the polarizability and the Green's tensor defined we now need to define the power spectrum and the quantum friction
```cpp
// define power spectrum
PowerSpectrumHarmOsc power_spectrum(&greens_tensor, &polarizability);

// define quantum friction
QuantumFriction friction(&greens_tensor, &polarizability, &power_spectrum);
```
Perfect!
We have defined every object there is to calculate quantum friction and have set all parameters.
All that is left to do is calculate the quantum friction and print out the result.
Let us therefore define the numerical parameters and
