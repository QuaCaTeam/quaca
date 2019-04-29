***********
Development
***********

The project is currently developed in a private git repository.
The master branch should at all times be functional, which is why currently merges shall only be made in the presence of all developers.
To gain access to the code, contact the maintainer.


Structure
=========

The project is organized as follows::

   QuaCa/
   ├── bin
   ├── data
   ├── doc
   ├── include
   ├── src
   ├── test
   ├── Makefile
   └── README.md

In *bin* we store the binaries.
In *data* lie the input as well as the output files.
In *doc* we store the ``.rst``-files to make the documentation as well as the documentation itself.
For more information see `Documentation`_.
In *src* are the source files, i.e. the ``.c`` as well as the ``.h`` files.
In *test* we store the tests.
Finally the ``Makefile`` controls the whole project and the ``README.md`` is for information when you see the git repository.

Testing
=======

For testing we use `Check <https://libcheck.github.io/check/>`_. 
A nice introduction can be found `here <http://www.ccs.neu.edu/home/skotthe/classes/cs5600/fall/2015/labs/intro-check.html>`_.

To test Quaca, just run::

   $ make test && ./bin/check_quaca

Each time you run ``check_quaca`` there is an extended test log generated, which is stored in ``test/test.log``.


Documentation
=============

You probably need to install sphinx and the read-the-docs-theme with::

   $ apt-get install python3-sphinx
   $ pip install sphinx_rtd_theme

The documentation is organized like this::

   doc/
   ├── _build
   │   ├── doctrees
   │   └── html
   ├── conf.py
   ├── development.rst
   ├── img
   │   ├── quaca_medium.png
   │   ├── quaca.png
   │   └── quaca_small.png
   ├── index.rst
   ├── installation.rst
   ├── Makefile
   └── usage.rst

The documentation is written in all the ``.rst``-files.
There is a useful `cheat sheet <https://github.com/ralsina/rst-cheatsheet/blob/master/rst-cheatsheet.rst>`_ for the syntax.
See also the `Sphinx documentation <http://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_.

Then you run::
   
   $ make doc

to build the documentation, which is then build into ``_build/html``.
The documentation is build as html by default, this can be changed in the Makefile.
