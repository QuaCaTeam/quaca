# Write your own main
In the first tutorial you calculated the non-contact friction of an atom above the surface of a Drude bulk material using an input file and the executable `Friction`.
QuaCa can, however, also be used like a library and you can write your own executables!
In this tutorial, we will redo the last tutorial, but this time we will implement everything in our own main file.
You will learn how to:

1. Register a new executable.
2. Write the main function.
3. Build and run your program.
4. Check the results.

## 1. Register a new executable
Let us have a look at the subfolders of `quaca/`:
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

You will see a folder called `app/` and a folder called `src/` amongst them.
The folder `src/` contains all the code for the Green's tensors, polarizabilities and all other quantities that are needed for the calculation of quantum friction.
This folder does, however, *not* contain a main function, which is needed to build an executable (i.e. a program we can actually execute).
These main functions are inside the `app/` directory.
Look inside the `app/` folder and you will find e.g. a subfolder called `Friction` and a file called `CMakeLists.txt`.

Each folder inside the `app/` folder should contain a `.cpp`-file containing representing main file. 
Hence, each folder corresponds to an executable program.
In the last tutorial you have already worked with the executable `Friction`. 
Let us now create a new one.

First create a new folder inside the `app/` folder and call it `Tutorial`.
Next go to the file `CMakeLists.txt` inside the `app/` directory and add the following line to the end of the file
```cmake
add_subdirectory("Tutorial")
```
This tells CMake to consider this folder when building the executables.
Next create a file inside the `Tutorial/` folder called `tutorial.cpp` containing the following lines of code:
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
Create a file called `CMakeLists.txt` inside the `Tutorial/` folder containing the following lines of code:
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

With this setup, we are ready to test whether our tutorial executable works.
Go to the `quaca/build/` directory and type into the command line
```bash
quaca/build> cmake ..
```
You have to do this again, because we created new files and CMake needs to know about them!
Then build the project by typing
```bash
quaca/build> make
```
You should now see your executable in the `quaca/bin/` directory, so let us execute it:
```bash
quaca/bin> ./Tutorial
The tutorial executable works!
quaca/bin>
```
Good job!
Let us now proceed by filling our main with more meaningful code.

## 2. Write the main function
First let us define the permittivity.
We can do so by typing
```cpp
#include "Quaca.h"
#include <iostream>

int main(int argc, char *argv[])
{
  // parameters for permittivity
  double omega_p = 9.0;
  double gamma = 0.035;

  // define permittivity
  auto permittivity = std::make_shared<PermittivityDrude>(omega_p, gamma);

  return 0;
};
```
You have now created an object of the type PermittivityDrude called `permittivity` with the appropriate parameters.
This object can be used for the calculation of, e.g., the numerical value of the permittivity and can also be put into other objects as we will now do.
Define the reflection coefficients and use the permittivity as an argument by adding to our code
```cpp
// define reflection coefficients
auto refl_coefficients = std::make_shared<ReflectionCoefficientsLocBulk>(permittivity);
```
We hereby defined reflection coefficients and told them to use the permittivity from above.
Let us proceed by defining our Green tensor using the reflection coefficients as an argument.
```cpp
// parameters for green's tensor
double v = 1e-4;
double beta = 1e6;
double z_a = 0.025;

// numerical error for green's tensor
double delta_cut = 20;
vec::fixed<2> rel_err = {1E-4, 1E-2};

// define the Green's tensor
auto greens_tensor = std::make_shared<GreensTensorPlate>(v, beta, z_a, refl_coefficients, delta_cut, rel_err);
```
Let us proceed by defining the polarizability.
```cpp
// parameters for polarizability
double omega_a = 1.3;
double alpha_zero = 6e-9;

// define polarizability
auto polarizability = std::make_shared<Polarizability>(omega_a, alpha_zero, greens_tensor);
```
With the polarizability and the Green's tensor specified we now need to define the power spectrum and the quantum friction
```cpp
// define power spectrum
auto power_spectrum = std::make_shared<PowerSpectrum>(greens_tensor, polarizability);

// numerical error for quantum friction
double rel_err_omega = 1e-1;

// define quantum friction
auto friction = std::make_shared<Friction>(greens_tensor, polarizability, power_spectrum, rel_err_omega);
```
Perfect!
We have now defined every object that is needed to compute quantum friction and have fixed all parameters.
All there is left to do is to calculate the friction force and print the result.

For the calculation we still need to define our scale.
We want to calculate the friction for velocities ranging from `1e-4` to `1e-2` on a logarithmic scale.
```cpp
// loop over v
double start = 1e-4;
double end = 1e-2;
int number_of_steps = 40;
double spacing = pow(end / start, 1. / ((double)number_of_steps - 1.0));
```
We also need to define an output file like this
```cpp
// define output file
std::ofstream file;
file.open("tutorial_mainfile.csv");
```

Now we put it all together in a loop
```cpp
double step, value;
for (int i = 0; i < number_of_steps; i++) {
  step = start * pow(spacing, i);
  friction.get_greens_tensor()->set_v(step);
  value = friction->calculate(NON_LTE_ONLY);

  file << step << "," << value << "\n";
};
```
The option `NON_LTE_ONLY` specifies the form of the power spectrum we want to use when calculating friction; for more information have a look at the [documentation](api/powerspectrum).

At the end of our main function we have to close the file, so that all together our `tutorial.cpp` looks like this
```cpp
#include "ProgressBar.hpp"
#include "Quaca.h"

int main(int argc, char *argv[]) {

  // parameters for permittivity
  double omega_p = 9.0;
  double gamma = 0.035;

  // define permittivity
  auto permittivity = std::make_shared<PermittivityDrude>(omega_p, gamma);

  // define reflection coefficients
  auto refl_coefficients = std::make_shared<ReflectionCoefficientsLocBulk>(permittivity);

  // parameters for green's tensor
  double v = 1e-4;
  double beta = 1e6;
  double z_a = 0.025;

  // numerical error for green's tensor
  double delta_cut = 20;
  vec::fixed<2> rel_err = {1E-4, 1E-2};

  // define the Green's tensor
  auto greens_tensor = std::make_shared<GreensTensorPlate>(v, beta, z_a, refl_coefficients, delta_cut, rel_err);

  // parameters for polarizability
  double omega_a = 1.3;
  double alpha_zero = 6e-9;

  // define polarizability
  auto polarizability = std::make_shared<Polarizability>(omega_a, alpha_zero, greens_tensor);

  // define power spectrum
  auto power_spectrum = std::make_shared<PowerSpectrum>(greens_tensor, polarizability);

  // numerical error for quantum friction
  double rel_err_omega = 1e-1;

  // define quantum friction
  auto friction = std::make_shared<Friction>(greens_tensor, polarizability, power_spectrum, rel_err_omega);

  // loop over v
  double start = 1e-4;
  double end = 1e-2;
  int number_of_steps = 40;
  double spacing = pow(end / start, 1. / ((double)number_of_steps - 1.0));

  // define output file
  std::ofstream file;
  file.open("tutorial_mainfile.csv");

  double step, value;
  for (int i = 0; i < number_of_steps; i++) {
    step = start * pow(spacing, i);
    friction->get_greens_tensor()->set_v(step);
    value = friction->calculate(NON_LTE_ONLY);

    file << step << "," << value << "\n";
  };

  // close file
  file.close();

  return 0;
};
```

## 3. Build and run your program
We are now ready to build and run what we have programmed.
Go into the `build/` directory and run
```bash
quaca/build> make Tutorial
```
You have now built the Tutorial executable which is automatically stored in the `quaca/bin/` folder.
Run the executable by typing into the console
```bash
quaca/build> ./../bin/Tutorial
```
Notice that, because of the way we have written our main function, the output file `tutorial_mainfile.csv` will be placed wherever the executable is run, so in this case it will be placed in the `build/` directory.

## 4. Check the results
Plot the data again by going to the `quaca/plots` directory and using our prepared plot script
```bash
quaca/plots> python plot.py ../build/tutorial_mainfile.csv
```
You should obtain a plot that looks exactly like the one we have gotten in the [first Tutorial](tutorials/first_calculation).
