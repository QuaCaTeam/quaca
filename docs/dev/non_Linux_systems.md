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
even if gcc/g++ is called, is under the hood the clang compiler. It used to be the case
that this version of the clang compiler does not support OpenMP natively. Eventhough,
a support can be achieved by installing the library *libomp*.
For simplicity, we will assume that [HomeBrew](brew.sh/index_de.html) is used for managing packages.
To install the new library simply type in the terminal
```asm
brew install libomp
```
This should install the library *libomp* in the directory: /usr/local/opt/libomp.
For some reason, at least in my configuration (macOS 10.14.6 Mojave) still CMake
is not able to detect the OpenMP support on its own. Therefore we have to pass
the correct directories to CMake, when we invoke it for the first time.
We start equivalent to Linux distributions by navigating in the terminal to the
*QuaCa* directory
```asm
cd path-to-your-local-quaca
```
Next we create the build directory and navigate to it.
```asm
mkdir build && cd build
```
Now we want to run cmake for the first time. Here we now have to pass all the releveant
variables to cmake, such that it detects the installed libomp library
```asm
cmake .. -DOpenMP_CXX_FLAGS="-preprocessor -fopenmp /usr/local/opt/libomp/lib/libomp.dylib -I/usr/local/opt/libomp/include -DOpenMP_CXX_LIB_NAMES="libomp"
```
Now a simple
```asm
make
```
should build the full project.
In principle we could have also changed the compiler used by cmake to a newer version of
clang which gets delivered by the llvm library or by installing the correct gcc/g++
compiler. Eventhough at least in the latter case we will run into the issue, that the installed
boost library was build with the default clang compiler which relies on a different
version of gcc then the newly installed gcc. Therefore, linking to boost will result
in compiler errors. 
The whole issue with OpenMP, Boost and macOS is nicely summarized in this [StackOverflow question](https://stackoverflow.com/questions/47081991/is-c-compilable-with-openmp-and-boost-on-macos/47225639#47225639)
For some reasons, the suggested fixes did not work for me. Therefore, I used the solution 
posted in this [StackOverflow question](https://stackoverflow.com/questions/60126203/how-do-i-get-cmake-to-find-openmp-c-openmp-cxx-etc).
For some reason in the given answers some flags are set several times which I omitted as well
as all options for the c-compiler, which are not needed for *QuaCa*.
