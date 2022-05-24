---
title: 'QuaCa: an open-source library for fast calculations of long-time atom-surface interactions'
tags:
  - C
  - C++
  - fluctuation-induced phenomena
  - atom-surface interaction
  - quantum friction
  - Casimir friction
  - Casimir effect
authors:
  - name: Marty Oelschläger
    affiliation: 1
  - name: Simon Hermann
    affiliation: 1
  - name: Christoph H. Egerland
    orcid: 0000-0002-1099-6433
    affiliation: 1
  - name: Daniel Reiche
    affiliation: "1, 2"
affiliations:
 - name: Humboldt-Universität zu Berlin, Institut für Physik, AG Theoretische Optik & Photonik, 12489 Berlin, Germany
   index: 1
 - name: Max-Born-Institut, 12489 Berlin, Germany
   index: 2
date: 4 June 2020
bibliography: paper.bib
---

# Summary

QuaCa is a package for the investigation of nonequilibrium long-time atom-surface interactions. One of its core applications is the computation of the noncontact friction of a microscopic object moving with constant velocity $v$ and height $z_a$ above a macroscopic surface along a uniform trajectory (see autoref{fig:setup}). QuaCa provides the direct calculation of this fully nonequilibrium force acting upon the microscopic object. This peculiar force was investigated and discussed by many theoretical publications [@dedkov2002; @pendry2010; @scheel2009; @volokitin2002], however, could not yet be experimentally confirmed. Recently, Farías et al. [@farias2020] suggested that the experimental measurement is on the verge of a potential breakthrough. Thus, the interest in precise and complete calculation of this force are twofold. On the one hand, further potential setups considering different materials and geometries can easily be explored in order to hint experimental feasible setups. On the other hand it can be used to evaluate potential future applications. 

![Sketch of the setup. A microscopic object, here depicted as an atom, moves with constant velocity and height abov a flat macroscopic surface..\label{fig:setup}](images/setup.svg)

The noncontact friction is caused by (quantum) fluctuations of the electromagnetic fields and is thus related to the Casimir-Polder or van der Waals force [@vanderwaals1873; @casimir1948]. In contrast to the Casimir-Polder force the noncontact friction is a pure nonequilibrium phenomenon, without an equilibrium counterpart [@intravaia2019]. Besides the complexity of a nonequilibrium phenomenon the frictional force exhibits a finite memory, called non-Markovian, with respect to the interaction with the fluctuating electromagnetic fields [@intravaia2016a].  Due to its finite memory, the evaluation of the friction involves a fivefold integration of non-trivial functions (including steep poles and oscillations of the integrands).

With QuaCa we implemented the efficient integration scheme outlined by Oelschläger [@oelschlager2019] and thus allow for a fast and efficient computation. By virtue of its modularity QuaCa can also be used as a library for several quantities used in the context of light-atom or light-matter interaction (as e.g. the power spectrum, the polarizability or the (integrated) Green's tensor). Therefore, other interesting physical mechanisms derived from these quantities, as e.g. the provided decay rate, can be easily calculated within the framework of relative motion and nonequilibrium.

The implementation of the noncontact friction evaluates the force in its nonequilibrium steady state. In the long-time regime, i.e. when all transients have decayed and the system reaches its nonequilibrium steady-state, the noncontact friction for a harmonic dipole moment can be solved completely without further perturbational approximations. For the zero temperature case Intravaia et al. [@intravaia2014;@intravaia2019] presented the compact expression

$$
  F_\mathrm{fric} = 
-2
\int_{0}^{\infty} \mathrm{d}\omega\, \int\frac{\mathrm{d}^2\mathbf{k}}{(2\pi)^2} \, k_x\,
\mathrm{Tr}\left[
\underline{S}^\intercal(\mathbf{k}^\intercal\mathbf{v}-\omega)\underline{G}_\Im(\mathbf{k}, z_a, z_a, \omega)
\right],
$$

where $\underline{S}(\omega)$ is the power spectrum of the dipole moment's autocorrelator, $\underline{G}(\mathbf{k},z_a,z_a,\omega)$ is the electric Green's tensor $\mathbf{k}=(k_x,\,k_y)^\intercal$ is the two-dimensional wavevector, and $\omega$ is the frequency. The power spectrum itself contains the polarizability $\underline{\alpha}(\omega)$ and a twofold integration over the Green's tensor and some weighting function. This framework was extended to rotational degrees of freedom [@intravaia2019] and to finite temperatures [@oelschlager2019] in recent years. These important extensions are implemented in QuaCa, which enables to calculate state-of-the-art results of this interesting nonequilibrium phenomenon.

From a coding perspective QuaCa features:

- A library implementing several important quantities used in light-matter interaction
- Two readily written executables:
  1) calculation of the noncontact friction
  2) calculation of the decay rate
- Use of modern integration routine cquad [@galassi2019] to handle the demanding integrals
- Support for parallelization with OpenMP
- Extensive unit and integration tests
- Use of modern C++

QuaCa supports following geometries, material models and memory kernels:
 - Semi-infinite bulk and finite slab geometry
 - Drude and Drude-Lorentz permittivity model
 - Ohmic and Single-Phonon (see [@lopez2018; @oelschlager2019 ]) memory kernel
 
When used as a library the following quantities are separately accessible:
 - the polarizability $\underline{\alpha}(\omega)$ (considering backaction of the fields onto the atom)
 - the (integrated) Green's tensors $\underline{G}(\mathbf{k},z,z'\to z,\omega)$ for the above described geometries
 - the power spectrum of the dipole moment's autocorrelator $\underline{S}(\omega)$
 - the reflection coefficients $r(\omega,\mathbf{k})$, permittivity $\epsilon(\omega)$, and memory kernel $\mu(\omega)$

Due to its modular nature QuaCa is easily extendable to more geometries, material models and memory kernels, which will be within the scope of future releases.

While the related Casimir force sparked excellent code packages [@hartmann2020;@reid2015;@oskooi2010] for its calculation, QuaCa is to the best of our knowledge the first package/library which provides an easy accessible framework to compute the noncontact frictional force.

QuaCa has already been used in several publication [@reiche2018;@intravaia2019;@oelschlager2018;@reiche2018] and provided novel scientific insight into the realm of nonequilibrium fluctuation-induced phenomena.

# Acknowledgements
The authors thank Bettina Beverungen, Kurt Busch, Dan-Nha Huynh and Francesco Intravaia for many fruitful discussions and an always open door.

# References