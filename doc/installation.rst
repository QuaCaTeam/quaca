************
Installation
************

To build QuaCa from source simply clone the git repository via::

    $ git clone git@git.physik.hu-berlin.de:egerlanc/quaca.git

Note that this repository is not public, so for now you need to ask the developers for access.
To build the project type::

    $ make quaca
    $ sudo make install

The default directory for installation is ``/usr/local/bin`` and can be changed in the *Makefile*.


Dependencies
============

GCC for compilation, GSL library for integration routines, Check for testing.

To install all dependencies run::

   $ sudo apt-get install gcc build-essential libgsl-dev check 

For the quantum friction in a cylinder we need to compute Bessel functions of a complex argument, which the GSL unfortunately does not provide.
We choose to use `Arb <http://arblib.org/>`_ for this.
To install Arb run::

   $ sudo apt-get install libflint-arb-dev
