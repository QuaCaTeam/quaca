# QuaCa <!-- {docsify-ignore-all} -->
**QuaCa** stands for Quantum / Casimir Friction and is a numerical framework to calculate just that.

## Dependencies
QuaCa requires:

* [GSL](https://www.gnu.org/software/gsl/): integration routines
* [Boost](https://www.boost.org/): file parser and file handling
* [BLAS](http://www.netlib.org/blas/) and [LAPACK](http://www.netlib.org/lapack/): fast linear algebra

Furthermore we use:
* [Armadillo](http://arma.sourceforge.net/): wrapper for linear algebra (uses LAPACK and BLAS)
* [Catch2](https://github.com/catchorg/Catch2): unit testing
* [docsify](https://docsify.js.org): documentation

Those prerequisites can easily be installed via `apt install` via the following command
```bash
sudo apt-get update && sudo apt-get install libblas-dev liblapack-dev cmake make g++ libboost-filesystem-dev libboost-program-options-dev
```

or with [conda](https://docs.conda.io/en/latest/):

```bash
conda env create --file environment.yml
conda activate quaca-env
```

We currently do not support installation on Mac or Windows.

## Installing
To obtain the source code type
```bash
git clone https://github.com/QuaCaTeam/quaca.git
```
into the console

Then to build the code type inside the project directory (if you use `conda` please activate `quaca-env`)
```bash
$ cd quaca
quaca $ mkdir -p build
quaca $ cd build
quaca/build $ cmake .. -DCMAKE_PREFIX_PATH=$CONDA_PREFIX #or just cmake .. if you don't use conda
quaca/build $ make

```

Two executables called `Friction` and `Decay` should now have been build and can be found in the `quaca/bin` directory.

## Documentation
For now we do not upload the documentation anywhere, so it can only be viewed locally using the following command from the command line
### Python 2.7
```bash
quaca $ cd docs && python -m SimpleHTTPServer 3000
```
### Python 3
```bash
quaca $ cd docs && python -m http.server 3000
```

You can then see the documentation in your browser at the address `http://localhost:3000/`.

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
More detailled information is given in the section [Testing](https://quacateam.github.io/quaca/#/dev/testing).

## Benchmark
For a benchmark we used a Intel(R) Core(TM) i7-7500U CPU @ 2.70GHz and Intel(R) Core(TM) i9-10900K CPU @ 3.70GHz and ran the tutorial configuration with 

```bash
quaca/bin $ ./Friction --file ../data/.tutorial.json --threads {threads}
```
and obtain following results

|Threads|max. 2.7GHz|max. 3.7GHz|
|---|---|---|
|  1|  1477 s | 623 s  |
|  2|  995 s |  319 s  |
|  3|  703 s |  222 s  |
|  4|  555 s  |  170 s |

## Contributing
Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on how you can become a part of this project.

## Authors
- [Marty Oelschläger](mailto:myoel@posteo.de)
- Simon Hermann
- Christoph Egerland
- Daniel Reiche
