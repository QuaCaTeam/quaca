# File input and API

## Command line options
QuaCa is envoked in the command line and uses [boost program options](https://www.boost.org/doc/libs/1_71_0/doc/html/program_options/tutorial.html) for its API.
The possible options are:

## Input File
Parameters are given to QuaCa in the form of `.ini` files.

An input file always follows the following format

```ini
[MemoryKernel]
type=ohmic
gamma=3
```

For a more hands on guide to the input files see the [Tutorial](tutorial).
