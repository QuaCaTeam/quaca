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

### Installing
To obtain the source code type into the console
```
git clone git@git.physik.hu-berlin.de:top/codeprojects/quaca.git
```

Then to build the code type inside the project directory
```
quaca $ mkdir -p build
quaca $ cd build
quaca/build $ cmake ..
quaca/build $ make quaca
```

An executable called `quaca` should now have been build and can be found in the `quaca/bin` directory.

### Documentation
The documentation can either be seen [here](https://QuaCaTeam.github.io/quaca) or be viewed locally using the following command from the command line
```
quaca $ cd docs && python -m SimpleHTTPServer 3000
```
You can then see the documentation in your browser at the adress `http://localhost:3000/`.

### Testing
The follow a test-driven development style and are currently use [Catch2] for our unit testing.
To run all unit tests, build the tests first with
```
quaca/build $ make test_quaca
```
then run the tests from the `bin/` directory
```
quaca/bin $ ./test_quaca
```

### Contributing
Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on how you can become a part of this project.

### Authors
- Marty Oelschläger
- Simon Hermann
- Christoph Egerland
