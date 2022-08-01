[![Build Status](https://travis-ci.com/myoelmy/quaca.svg?branch=master)](https://travis-ci.com/myoelmy/quaca)

# QuaCa
**QuaCa** stands for Quantum / Casimir Friction and is a numerical framework to calculate just that.

## Getting Started

### Dependencies
QuaCa requires for its installation:

* [GSL](https://www.gnu.org/software/gsl/): integration routines
* [Boost](https://www.boost.org/): file parser and file handling
* [BLAS](http://www.netlib.org/blas/) and [LAPACK](http://www.netlib.org/lapack/): fast linear algebra

Furthermore we use:
* [Armadillo](http://arma.sourceforge.net/): wrapper for linear algebra (uses LAPACK and BLAS)
* [Catch2](https://github.com/catchorg/Catch2): unit testing
* [docsify](https://docsify.js.org): documentation

## Installing
To obtain the source code type
```bash
git clone https://github.com/QuaCaTeam/quaca.git
```
into the console.

Then to build the code type inside the project directory (if you use `conda` please activate `quaca-env`)
```bash
$ cd quaca
quaca $ mkdir -p build
quaca $ cd build
quaca/build $ cmake .. -DCMAKE_PREFIX_PATH=$CONDA_PREFIX #or just cmake .. if you don't use conda
quaca/build $ make

```

Two executables called `Friction` and `Decay` should now have been build and can be found in the `quaca/bin` directory.

### Documentation
The documentation can either be seen [here](https://quacateam.github.io/quaca/) or be viewed locally using [Python 3](https://www.python.org/download/releases/3.0/) by
```
quaca $ cd docs && python3 -m http.server 3000
```
or [Python 2.7](https://www.python.org/download/releases/2.7/) by
```
quaca $ cd docs && python -m SimpleHTTPServer 3000
```
You can then see the documentation in your browser at the adress `http://localhost:3000/`.

## Testing
Employing test-driven development, we are using [Catch2](https://github.com/catchorg/Catch2) for our unit and integrated testing.
To build and run all implemented test use
```bash
quaca/build $ make test_quaca_unit
quaca/build $ make test_quaca_integrated
```
and afterwards run the tests from the `bin/` directory
```bash
quaca/bin $ ./test_quaca_unit
quaca/bin $ ./test_quaca_integrated
```
More detailled information is given in the section [Testing](dev/testing.md).

### Contributing
Please read [contributing](CONTRIBUTING.md) for details on how you can become a part of this project.


### Authors
- Marty Oelschläger
- Simon Hermann
- Christoph Egerland
- Daniel Reiche
