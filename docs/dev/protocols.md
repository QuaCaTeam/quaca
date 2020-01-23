Protocols
-----

### Meeting 1: Christoph, Simon and Marty (14.01.2020)

- Each participant presented his progress
    * Implementation of integration options as demanded accuracy were discussed
    * Relative error bounds are given as input parameter in the .ini file
    * Each Class draws its respective error bound from the .ini file

- Testing was discussed
    * Marty will investigate with respect to physical domains of each input parameter
    * The catch2 approximate comparison shall be done by Approx(...).epsilon(...) to investigate relative errors (Approx(...).margin(...) sets absolute error bounds, which might not be reasonable)

- Branch structure was discussed
    * Christoph proposed to set up a __developement__ branch
    * Currently the __master__ branch represents a __developement__ branch
    * When the first stable version is released, we should switch to a __master__+__developement__ branch structure

- As an outlook following task were assigned:
    * Christoph continues to test the polarizability
    * Simon continues to test the free space friction and the power spectrum
    * Marty's tasks are ...
        * rewrite the Relations PDF
        * add the &omega;<sub>a</sub> << &omega;<sub>cut</sub> case for an object with internal bath
        * test the GreensPlate class
        * investigate physical domains of each input parameter