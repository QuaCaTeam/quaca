---
title: 'QuaCa: an open-source library for fast calculations of steady-state quantum friction'
tags:
  - C
  - C++
  - fluctuation-induced phenomena
  - atom-surface interaction
  - quantum friction
  - Casimir friction
authors:
  - name: Marty Oelschläger
    affiliation: "1,2,3"
    orcid: 0000-0002-0120-9733
  - name: Simon Hermann
    affiliation: 3
  - name: Christoph H. Egerland
    orcid: 0000-0002-1099-6433
    affiliation: 3
  - name: Daniel Reiche
    affiliation: 3
    orcid: 0000-0002-6788-9794
affiliations:
 - name: Max-Born-Institut, 12489 Berlin, Germany
   index: 1
 - name: dida Datenschmiede GmbH, Hauptstraße 8, 10827 Berlin, Germany
   index: 2
 - name: Humboldt-Universität zu Berlin, Institut für Physik, 12489 Berlin, Germany
   index: 3
date: 25 May 2023
bibliography: paper.bib
---

# Summary

QuaCa is an extensible library facilitating the computation of steady-state atom-surface quantum friction.
Due to its modular domain-driven structure, QuaCa can be of further use to calculate relevant quantities that are often needed in the context of electromagnetic dispersion forces.
Quantum (or Casimir) friction is a quantum-optical fluctuation-induced force that occurs in dynamical nonequilibrium, i.e. when a number of bodies are moving relative to one another [@pendry1997;@scheel2009].
The frictional interaction between the moving system of interest and its environment is mediated by the (material-modified) quantum vacuum and persists even at zero temperature.
At finite temperatures, thermal fluctuations can modify the interaction and quantum friction can be connected to the Einstein-Hopf effect [@einstein1910;@milonni1981;@ford1985;@mkrtchian2003;@oelschlaeger2021].
For a comprehensive review of quantum friction in various contexts as well as related effects in dynamical nonequilibrium, we refer, for example, to the reviews by @volokitin2007, @dedkov2017, @reiche2022 and the references therein.

# Statement of need

Due to the relative weakness of quantum friction with respect to related fluctuation-induced effects in equilibrium, such as the Casimir(-Polder) or the van der Waals effect [@dalvit2011], quantum friction has yet evaded experimental confirmation.
This has instigated a surge of interest in exploring potentially useful experimental designs [@volokitin2011;@milton2016;@farias2020;@lombardo2021;@reiche2022]. To promote the effect to the measurable realm a sophisticated (numerical) optimization of the experimental setup appears necessary [@oelschlager2019;@reiche2021]. 
This need can be met using QuaCa. 
QuaCa, to the best of our knowledge, pioneers the numerical simulation of quantum friction and, to date, is the only ready-to-use and openly available package for exploring experimental designs that can be utilized to measure quantum friction.


# Physics of quantum friction

We focus on the nonequilibrium steady-state of a microscopic particle moving with non-relativistic velocity $v$.
The frictional force connected to the electric interaction then reads [@intravaia2014;@intravaia2016;@intravaia2016a;@intravaia2019]
$$
  F_\mathrm{fric} =
-2
\int_{0}^{\infty} \mathrm{d}\omega\, \int\frac{\mathrm{d}h}{2\pi} \, h\,
\mathrm{Tr}\left[
\underline{S}^\intercal(hv-\omega)\underline{G}_\Im(h, \mathbf{R}_a, \omega)
\right],
$$
where $\underline{S}(\omega)$ is the particle's power spectrum, $\underline{G}(h,\mathbf{R}_a,\omega)$ is the electric Green tensor, $h$ is the wavevector along the direction of motion, $\mathbf{R}_a$ gives the position of the particle in the plane perpendicular to the direction of motion, $\omega$ is the frequency, the superscript $\intercal$ denotes the transpose of a matrix and $\underline{G}_{\Im}=(\underline{G}-\underline{G}^{\dagger})/(2i)$.
The power spectrum encodes the temperature-dependence of the force as well as the polarizability $\underline{\alpha}(\omega)$ of the microscopic particle.
The latter is dressed by the electromagnetic environment -- again given by the electric Green tensor.
This approach respects backaction from the environment onto the system [@intravaia2016;@reiche2020a] and includes finite temperatures [@oelschlager2019;@oelschlaeger2021] as well as any net transfer of spin angular momentum [@intravaia2019].


![Schematic of the physical model. A microscopic object, for example, an atom, moves with constant velocity and height above a flat macroscopic surface. 
The particle is attracted by the surface ($F_{\mathrm{CP}}$, Casimir-Polder force) and experiences a moderating force ($F_{\mathrm{fric}}$, quantum friction). \label{fig:setup}](images/setup.svg){ width="800" height="600" style="display: block; margin: 0 auto" }

# Numerical approach and code structure

QuaCa (i) allows for a computation of $F_{\mathrm{fric}}$ on any regular personal computer and (ii) facilitates extensibility due to its modular code structure.

(i) Computing quantum friction in the form presented above, the biggest challenge arises from the evaluation of nested frequency- and wavevector integrals.
Further, both the polarizability and the Green tensor feature a number of poles that can lead to fastly oscillating integrands.
The application of a Wick rotation that transforms the oscillations into exponential decays, as it became common in the numerical treatment of equilibrium fluctuation-induced effects [@oskooi2010;@johnson2011;@reid2015;@hartmann2020], is non-viable due to the Doppler-shift of the frequency.
Here, QuaCa uses an analytically equivalent version of $F_{\mathrm{fric}}$, which is particularly suited for numerical computations. We decouple nested integrals and take explicit care of the occurring poles. For details of the procedure, we refer to @oelschlager2019.

(ii) We have chosen a modular implementation strategy guided by physical principles. 
We assigned each meaningful physical observable in the setup to individual objects in the package which can be tested, replaced, computed or used independently. 
This allows for a simple adaption to different geometries or materials and physically relevant objects of the code can be used as a library.
In the current version, we demonstrate this versatility by including a routine for computing the atomic decay rate.
Numerical optimizations in the context of design and inverse-design [@molesky2018;@bennett2020] using appropriate Maxwell solvers [@busch2011] as an input to QuaCa can be one future application of the modular code package.


QuaCa contains the following main features and characteristics:

- A library implementing several important quantities used in light-matter interaction
- Two ready-to-use executables:
  1) calculation of the friction force
  2) calculation of the decay rate
- Use of the integration routine cquad [@galassi2019] 
- Support for parallelization with OpenMP
- Extensive unit and integration tests
- Use of modern C++

QuaCa currently supports the following geometries, and material models:

 - Semi-infinite bulk and finite slab geometry
 - Drude and Drude-Lorentz permittivity model
 - Ohmic and single-phonon [@lopez2018;@oelschlager2019] memory kernel for possible internal degrees of freedom of the moving particle

When used as a library the following quantities are separately accessible:

 - the dressed polarizability
 - the Green tensors for the above described geometries
 - the power spectrum of the particle's dipole moment
 - the reflection coefficients of a planar interface and the permittivity

QuaCa has been used to study the impact of multilayer structures [@oelschlager2018] or spatial non-locality in the bulk material [@reiche2019] on quantum friction, to investigate the net angular momentum transfer between a moving atom and its electromagnetic environment [@intravaia2019] and to explore the role of finite temperatures in the context of the thermal viscosity of the material-modified vacuum [@oelschlaeger2021]. The package is released under the MIT-license.

# Acknowledgements
The authors thank Bettina Beverungen and Dan-Nha Huynh for fruitful discussions and are particularly grateful to Kurt Busch and Francesco Intravaia for an inspiring guidance during all stages of the project.

# References
