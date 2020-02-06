# Contributing

Thank you for your interest in contributing to this project!
Below you'll find all information regarding the development of QuaCa.

We follow an agile-ish development style utilizing test-driven development and continuous integration.
All changes should be first discussed with the maintainers either via email or in a personal conversation.

## Workflow
The project features a master branch which always contains a working version of QuaCa.
New features are implemented in feature branches, which, if all tests succeed, are merged into the master.
This also outlines the development style for any new developer: branch out the master, implement your stuff, test it, make adjustments and after some communication with the maintainers merge it into the master.

## Documentation
There are two parts in which QuaCa is documented.
The first part are the comments within the code, explaining how exactly certain things are achieved.
The second part consists of a user friendly documentation in tutorial style and a developers documentation giving more insight into how QuaCa works under the hood.

## Coding conventions and styling
We follow the [C++ Core Guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines).
Some important conventions we follow are:
* Use a `.cpp` suffix for code files and `.h` for interface files.
* Use include guards for all `.h` files.
* Use a separate `.h` file for each class.
* Class names are written in capital style (`House`, `GreensTensor`, `AnEvenLongerName`, ...)
* Variable and function names are written in underscore style (`calculate_advanced`, `get_name_of`, ...)
* Pointers should be abbreviated with `pt`

For styling we use [clang-format](https://clang.llvm.org/docs/ClangFormat.html).
Also have a look [here](https://electronjs.org/docs/development/clang-format) on how to implement clang-format in your editor.
