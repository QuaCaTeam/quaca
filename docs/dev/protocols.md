Protocols
-----

### Meeting 1: Christoph, Simon and Marty (14.01.2020)

- Each participant presented his progress
    * Implementation of integration options as demanded accuracy were discussed
    * Relative error bounds are given as input parameter in the .json file
    * Each Class draws its respective error bound from the .json file

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

### Meeting 2: Christoph, Simon and Marty (21.01.2020)

- Christoph introduced clang format

- Simon reported on testing issues (power spectrum is too much zero, and random() produces always the same values)

- Marty presented code progress on GreensTensorPlate

- Notes one testing relations have been discussed. Improvement is needed

- Next week's goal: Calculate Friction

- Christoph presented Travis for future move of the project

- Discussed the handing over to Francescos:
    * Francesco will be responsible for the survival of the code
    * It is a very easy possibility to produce publications and generate projects for students
    * Christoph, Simon and Marty continue to assure the quality of the code for a transition period

### Meeting 3: Marty and Christoph (04.02.2020)

- all showed their code progress
- we need a class for summing Green's tensors
- quantumfriction should be an abstract class
- marty keeps working on the reflection coefficients
- quantum friction with plate setup works and will soon be pushed into the master
- optimization and code profiling should be made a standing issue

### Meeting 4: Christoph, Simon and Marty (09.04.2020)
 - Stand der Dinge + Kapazität
    * Simon uebergibt lead development an bettina
    * kapazitaet ist da

 - Was wollen wir noch umsetzen?
    * Umstrukturierung der Test Cases (kein random testing, lieber konkrete parameter)
    * Mehr Groessen bestimmen (Power spectrum, Drehmoment, ...)
    * Dokumentation polieren


  - Wo veroeffentlichen wir?
    * Marty fragt Kurt bzgl JOSS


  - Was ist fuer eine Veroeffentlichung zu tun?
   * Code auf guten Stand bringen
   * Porten auf GitHub


  - Anderes Paper
   * Konkrete Aufgaben kommunizieren
   * Marty und Daniel kümmern sich um Inhalte
   * Christoph und Simon helfen beim debugging von QuaCa


  - Vorgehensweise
    * Wöchentliche "Treffen" (Donnerstag 10 Uhr)
    * _Marty_: Dokumentation, Kommunikation mit Kurt
    * _Simon_: Refactoring von Test cases
    * _Christoph_: CMake und Memory Check, Profiling, Coverage

 - Timeline
   * Ende April: Code auf GitLab fertig
   * Ende Mai: fertig
   * 
   



Quaca Meeting 16.04.2020

1) Current state 
Christoph

    - from old to new school C++ (introduced shared_ptr in order to stop memory leaks)

    - openmp -> mpi


Marty (Documentation)

    - we now have a cover page

    - sections are clickable (link to first subsection)

    - sections and subsections were renamed

    - overall clean up

    - for permittivity and polarizability: classes with bath and nobath shall be unified by overloading the constructor

    - implementation of a unit conversion script

    - list of test cases were added to the documentation


Simon (Test cases)

    - for each relevant variable one test in each relevant magnitude


2) Paper & Choice of Journal

    - Marty prepares a mail for Kurt


3) Integration

    - Christoph is looking into another solution for the numerical integration (using templates)


4) App vs. Library

    - Options class will be moved to the executable


5) weitere executables?

    - would be nice

    - Marty: investigate potential executables (Casimir Polder, Angular Momentum, Free Energy ?)


-) TODO: 
Marty: - more Documentation
       - find executables
       - mail to Kurt
       - log lin ?
Christoph: - print_info() function

     - MPI

           - (delete options class and implement logger)
Simon: Test cases: rewrite GENERATE and pointer to smart pointer
        Example input files, Tutorials, print_info() function with comments in .csv file, Logger
        
        
# Quaca Meeting 24.04.2020

## Update

 - Simon Testcases: Test if tensors are different from zero

 - is_zero() function of the armadillo library is not working

    * We need to update armadillo library (to >= 9.850)


 - Christoph explains pointer problems in integrations:

    * Shared pointer cause massive overhead (with the current integration solution)

    * Integration could be adapted, but seems to be challenging


   - Simon is in favour of reimplementing, since it would improve the code's extensibility

 - Casimir-Polder and other features:

    * either we should spend some time on Casimir-Polder or leave it out

    * also Bettina plans to work on Casimir-Polder either way as a starting point of here PhD Thesis

    * alternatively we can implement the decay rate, well established quantity in the community and easy to implement

    (everything is already there)

    * rather implement angular momentum (good for presenting numerical efficiency of the integrations)


## Perspectives

    - Currently Bettina will not join the project

    - Paula want to do her thesis with QuaCa (September - December)

    - Daniel wants to learn the code

    - Is Francesco ready to hand over the code to another student?

    - Huddle with Francesco?


## Open Source Paper 

    - Kurt agrees to publish in JOSS

    - the developers are happy


## Tasks

    - __Current__:

     - Marty:

      - finish the docs

      - implement decay rate


     - Simon: 

      - finish tests


     - Christoph:

      - move options struct

      - openmp to mpi


    - __Coming up next__:


     - Marty:

      - Test private gitHub with Travis CI

      - Mail to Francesco


     - Simon: 

      - finish tests


     - Christoph:

      - shared pointer with new integration
