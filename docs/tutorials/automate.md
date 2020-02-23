!> TODO: Finish writing.

# Automate your calculation
As you saw in the last tutorial it is quite tedious to define all parameters in the main function itself.
It is also impractical, because if you want to change just one parameter, you have to recompile the whole project.
That's why in this tutorial we teach you how to automate these calculation, which means we show you how to use the `.json` files you have seen in the first tutorial.

In this tutorial we want to write an executable that prints the x-x-component of the polarizability tensor into a file for varying values of the frequency omega.
For this we have to do the following steps:

1. Register a new executable.
2. Write the main function.
3. Write a `.json` file.
4. Run the calculation.

## 1. Register a new executable
Create a folder called `AutoPol` inside the `quaca/app/` directory.
Add the following line to the end of the file `CMakeLists.txt` inside the `app/` directory
```cmake
add_subdirectory("AutoPol")
```
Create a file called `autopol.cpp` inside the `AutoPol/` folder and type inside it
```cpp
#include "Quaca.h"
#include <iostream>

int main(int argc, char *argv[])
{
  std::cout << "The autopol executable works." << std::endl;
  return 0;
};
```
Now create a file called `CMakeLists.txt` inside the `app/AutoPol/` folder and write
```cmake
# add executable
add_executable(AutoPol
  ${CMAKE_CURRENT_SOURCE_DIR}/autopol.cpp
  )
target_include_directories(AutoPol PRIVATE ../include)

# link libraries
target_link_libraries(AutoPol
  quaca
  )
```
Test the executable by building the project and then running the program
```bash
quaca/build> cmake .. && make
quaca/build> cd ../bin
quaca/bin> ./AutoPol
The autopol executable works.
quaca/bin>
```
Good job!
In the next step we will fill our main function.

## 2. Write the main function
We want to construct all of our classes using an input file.
QuaCa offers a class called `Options` which reads the input file given via the command line.
Write into your main function in `AutoPol/autopol.cpp`
```cpp
#include "Quaca.h"
#include <iostream>

int main(int argc, char *argv[])
{
  // get file from command line
  Options opts(argc, argv);
  std::string input_file = opts.get_parameter_file();

  return 0;
};
```

Then you can define all needed quantities, which in this case is only the polarizability with the input file.
The main function then looks like this.
```cpp
#include "Quaca.h"
#include <iostream>

int main(int argc, char *argv[])
{
  // get file from command line
  Options opts(argc, argv);
  std::string input_file = opts.get_parameter_file();

  // define polarizability
  PolarizabilityBath polarizability(input_file);

  return 0;
};
```
