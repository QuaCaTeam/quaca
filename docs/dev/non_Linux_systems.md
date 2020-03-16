# Using quaca in non Linux operating systems

This code project was mainly developed on Linux based operating systems.
Therefore, compiling and building the project in other build environments can be challenging.
In this section we want to collect functioning workarounds to get quaca
running also on other operating systems.

## macOS

In this section we will assume that all necessary library such as [boost](boost.org),
[GSL](gnu.org/software/gsl) and [LAPACK](netlib.org/lapack/) or [BLAS](netlib.org/blas/) have been installed.
While the prominent Linux distributions such as Ubuntu use the GNU gcc compiler as a standard compiler
this is not the case for macOS. The compilier installed via the XCode command line tool is,
even if gcc/g++ is called, under the hood the clang compiler. Specifically it is a version of
clang which does not support OpenMP. Therefore, naively trying to compile the *QuaCa* project with
CMake will not work.
In principle, there seem to be two ways out. Either one can install the proper gcc compiler
or one can install a clang compiler with OpenMP support. Generally, installing any other compiler with
OpenMP support should be an option as well.
For simplicity, we will assume that [HomeBrew](brew.sh/index_de.html) is used for managing packages.
To install either of the above compilers, simply run the command
```asm
brew install gcc
```
or
``` asm
brew install llvm
```
.
Next we have to go to the *QuaCa* directory type the following commands to create a directorz
for the build files and to change to that directory
```asm
mkdir build && cd build
```
If such a directory already exist it is better to delete it first, as cmake does not seem
to adapt to changes after it was run for the first time.
Now, beeing in the new build directory we have to initialize cmake. As cmake will use the default
c++ compiler of the operating system, we have to tell it, that we want to use a different one.
This can be achieved by the flag D followed with an option for the c++ compiler. The command should
look like this
```asm
cmake -D CMAKE_CXX_COMPILER=usr/local/bin/g++-9 ..
```
Here we assumed, that the version installed by HomeBrew was version 9. Please check which version
HomeBrew installed via
```asm
which g++-
```
If you want to use the new Clang compiler the cmake commakd should look like this
```asm
cmake -D CMAKE_CXX_COMPILER=usr/local/opt/llvm/bin/clang++ ..   
```
For some reason, on my system (macOS 10.14.6 Mojave) the first version with a new gcc compiler did not work
as there where some compiler errors raised by the boost library. I am not sure if this is connected
with an outdated version of the boost library or what else could have raised this issue. Eventhough,
the second version with an updated clang compiler works fine.
After this initial set-up there is no additional work to do, as cmake stores the path to the new compiler.
Only if the build directory is deleted for some reason the initial cmake commad as to follow the above syntax.
Furthermore, in this section I assumed that the new compilers where stored in the default directories.
Else the path to the directories have to be altered accordingly.

