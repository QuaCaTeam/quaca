*************
Documentation
*************

Usage
=====

To use QuaCa type into the command line::

    $ quaca INPUTFILE

So if you have an input file called ``params.in`` in a folder called ``data`` you would for example type into the terminal::

    $ quaca data/params.in

QuaCa then calculates its stuff and puts the values into an according .out file.

A sample .in file is depicted here

.. code-block:: bash

    # atomic constants
    v       5E-4   # in c
    za      5E-9   # in nm
    muquest 1      # if 1, the internal bath is included
    gamMu   1E-1   # in eV
    wa      1.3E0  # in eV
    a0      6E-9

    # material parameters
    einf  1.
    wp1   9.       # in eV
    g1    0.035E0  # in eV
    T     1E-6      # in K

    # numerical specifications
    kcut    30.
    relerr  1E-2
    recerr  1E-2
    abserr  1e-200

    # Range
    runvar v
    start  1E-4
    stop   1E-2
    steps  5
    scale  log



Analytics
=========

.. function:: double F0(void *p)

    Given a parameter struct this function computes the normalization constant

    .. math::

        F_0 = -3 \frac{\omega_{\mathrm{sp}}^5 \alpha_0}{2 \pi \varepsilon_0}


.. function:: double Fanat(void *p)

    Given a parameter struct this function computes the analytic limit of the translational contribution to quantum friction for small velocities of the atom.
    Notice that there are different results for different :math:`\mu`!

.. function:: double Fanar(void *p)

    Given a parameter struct this function computes the analytic limit of the rolling contribution to quantum friction for small velocities of the atom.


.. function:: double Ffreet(void *p)


.. function:: double Ffreer(void *p)

