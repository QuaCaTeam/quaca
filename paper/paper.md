---
title: 'QuaCa: an open-source library for fast calculations of long-time atom-surface interactions'
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
date: 23 June 2022
bibliography: paper.bib
---

# Summary

QuaCa is an extensible library facilitating fast computation of steady-state atom-surface quantum friction.
Quantum (or Casimir) friction is a quantum-optical fluctuation-induced force that occurs in dynamical nonequilibrium, i.e. when a number of bodies are moving relatively to one another [@pendry1997;@scheel2009].
The frictional interaction between the moving system of interest and its environment is mediated by the (material-modified) quantum vacuum and persists even at zero temperature.
At finite temperatures, thermal fluctuations can modify the interaction and quantum friction can be connected to the Einstein-Hopf effect [@einstein1910;@milonni1981;@ford1985;@mkrtchian2003;@oelschläger2021].
For a comprehensive review of quantum friction in various contexts as well as related effects in dynamical nonequilibrium, we refer, for example, to the reviews in Refs. [@volokitin2007;@dedkov2017;@reiche2022] and the references therein.

Due to the relative weakness of quantum friction with respect to comparable equilibrium effects, such as the Casimir(-Polder) or the van der Waals effect [@dalvit2011], quantum friction has, to the best of our knowledge, yet evaded experimental confirmation.
This has instigated a surge of interest in exploring potentially useful designs [@volokitin2011;@milton2016;@farias2020;@lombardo2021] (see also [@reiche2022] and references therein) and a sophisticated (numerical) optimization of the experimental setup appears necessary in order to promote the effect to the measurable realm [@oelschläger2019;@reiche2021].

That is where QuaCa can come into use.
QuaCa computes the quantum frictional force experienced by a microscopic particle (e.g. atom or nano-particle) moving along the invariant direction of an ensemble of macroscopic bodies.
We focus on the late-time regime where all the transitional dynamics has settled and the particle can be assumed to be in a nonequilibrium steady-state moving with constant non-relativistic velocity $v$.
In the regime of linear coupling between system and environment and in accordance with Refs. Intravaia et al. [@intravaia2014;@intravaia2016;@intravaia2016a;@intravaia2019], the frictional force connected to the electric interaction reads
$$
  F_\mathrm{fric} =
-2
\int_{0}^{\infty} \mathrm{d}\omega\, \int\frac{\mathrm{d}h}{2\pi)} \, h\,
\mathrm{Tr}\left[
\underline{S}^\intercal(hv-\omega)\underline{G}_\Im(h, \mathbf{R}_a, \omega)
\right],
$$
where $\underline{S}(\omega)$ is the power spectrum of the particle's dipole moment, $\underline{G}(h,\mathbf{R}_a,\omega)$ is the electric Green tensor specifying the (material-)properties of the environment, $h$ is the wavevector along the direction of motion, $\mathbf{R}_a$ gives the position of the particle in the plane perpendicular to the direction of motion (usually connected to the particle-surface distance) and $\omega$ is the frequency.
The power spectrum encodes the temperature-dependence of the force as well as the polarizability $\underline{\alpha}(\omega)$ of the microscopic particle.
The latter, due to the self-consistency of our approach [@intravaia2016;@reiche2020a], is dressed by the electromagnetic environment -- again given by the electric Green tensor.
This approach is fully respects backaction from the environment onto the system (exact in all orders coupling) and includes spin-momentum locking of confined light, which can lead, e.g., to a net transfer of spin angular momentum [@intravaia2019], as well as finite temperatures [@oelschlager2019;@oelschlaeger2021].

![Sketch of the setup. A microscopic object, say, an atom, moves with constant velocity and height above a flat macroscopic surface. At constant velocity $v$, the particle is attracted by the surface ($F_{\rm CP}$, Casimir-Polder force) and experiences a moderating force ($F_{\rm fric}$, quantum friction).. \label{fig:setup}](images/setup.svg)


The QuaCa package (i) allows for a fast computation of $F_{\mathrm{fric}}$ on any regular personal computer and (ii) facilitates efficient extensibility due to its modular code structure.

On the one hand (i), computing the quantum friction force in the form presented above, the biggest challenge arises from the evaluation of nested frequency- and wavevector integrals.
The nesting arises from the self-consistency of the approach and the Doppler-shift of the radiation.
Further, since our formalism takes realistic (dispersive and dissipative) materials into account, both the polarizability and the Green tensor feature a number of (physical) poles that can lead to fastly oscillating integrands, especially in the regime where retardation due to the finite speed of light comes into play.
The application of a Wick rotation that transforms the oscillations into exponential decays, as it became common in the numerical treatment of equilibrium fluctuation-induced effects [@oskooi2010;@johnson2011;@reid2015;@hartmann2020], is non-viable due to the Doppler-shift of the frequency.
Here, QuaCa uses an numerically optimized version of $F_{\mathrm{fric}}$, where some of the nested integrals decouple, puts particular care on the occurring poles, and hence allows for a fast and efficient computation of the force.
For details of the procedure, we refer to Ref. [@oelschlaeger2019].

On the other hand (ii), we have chosen an implementation strategy where any of the physical quantities, such as the permittivity, the polarizability, the Green tensor or the power spectrum, are assigned to individual objects in the package. These can be tested, replaced, computed or used independently.
By virtue of this modular structure, the QuaCa package can be easily adapted to different geometries and materials, parts of the code can be used as a library performing a subroutine in a larger project, or it can even be extended to compute other observables of fluctuation-induced light-matter interactions.
In the current version, we intended to demonstrate this versatility by including a routine for computing the atomic decay rate.
Numerical optimizations in the context of design and inverse-design [@molesky18;@bennett2020] using appropriate Maxwell solvers [@busch2011] can be one future application of the modular code package.


QuaCa contains the following main features and characteristics:

- A library implementing several important quantities used in light-matter interaction
- Two ready-to-use executables:
  1) calculation of the noncontact friction
  2) calculation of the decay rate
- Use of the integration routine cquad [@galassi2019] to handle the integrals
- Support for parallelization with OpenMP
- Extensive unit and integration tests
- Use of modern C++

QuaCa currently supports the following geometries, and material models:
 - Semi-infinite bulk and finite slab geometry
 - Drude and Drude-Lorentz permittivity model
 - Ohmic and single-phonon (see [@lopez2018;@oelschlager2019]) memory kernel for possible internal degrees of freedom of the moving particle

When used as a library the following quantities are separately accessible:
 - the dressed polarizability
 - the Green's tensors for the above described geometries
 - the power spectrum of the particle's dipole moment
 - the reflection coefficients of a planar interface, permittivity, and memory kernel for possible internal degrees of freedom

QuaCa has been used to study the impact of multilayer structures [@oelschlaeger2018] or spatial non-locality in the bulk material [@reiche2019] on quantum friction, to investigate the net angular momentum transfer between a moving atom and its electromagnetic environment [@intravaia2019] and to explore the role of finite temperatures in the context of the thermal viscosity of the material-modified vacuum [@oelschlaeger2021]. The package is released under the MIT-license.

# Acknowledgements
The authors thank Bettina Beverungen and Dan-Nha Huynh for fruitful discussions and are particularly grateful to Kurt Busch and Francesco Intravaia for an inspiring guidance during all stages of the project.

# References
