# Testing {docsify-ignore-all}

In order to ensure a solid quality of QuaCa we thoroughly implemented unit test for each building block, as well as integrated test. Besides usual (but still important) constructor tests which check that every object is transferred as intended, we tested following physical properties of the respective functions. Thus, we urge new contributor to write similar test suites for each new feature. 

## Green's Tensor
 - Reciprocity (for detailed information see for example [this wiki article](https://en.wikipedia.org/wiki/Reciprocity_(electromagnetism))) 
   * $\underline{G}(\omega,-\mathbf{k},z) = \underline{G}^\intercal(\omega,\mathbf{k},z)$
 - Reality of the quantity in time yields the crossing relation in frequency space
   * $\underline{G}(-\omega,\mathbf{k},z) = \underline{G}^\dagger(\omega,\mathbf{k},z) $
 - Asymptotic behavior for known regimes
   * For the [vacuum Green's tensor](api/greenstensor?id=greenstensorvacuum) in the limit of $v\ll c$ we find
    - $\int\frac{\mathrm{d}^2\mathbf{ k}}{(2\pi)^2}  \mathrm{Im}\left\{\underline{G}_0(\mathbf{ k}, z\to z',\omega+\mathbf{k}^\intercal\mathbf{v})\right\} \sim   \mathrm{diag}[1,1,1]\frac{\omega^3}{6\pi \epsilon_0 c^3}$

    - $\int\frac{\mathrm{d}^2\mathbf{ k}}{(2\pi)^2}  \mathrm{Im}\left\{\underline{G}_0(\mathbf{ k}, z\to z',\omega+\mathbf{k}^\intercal\mathbf{v})\right\}k_x \sim \frac{v/c }{6\pi\epsilon_0} 
 \frac{ \omega^4 }{c^4}   
\mathrm{diag}\left[
    1,\,
    2
    ,\,
     2\right]$

    - $\int\frac{\mathrm{d}^2\mathbf{ k}}{(2\pi)^2}  \mathrm{Im}\left\{\underline{G}_0(\mathbf{ k}, z\to z',\omega+\mathbf{k}^\intercal\mathbf{v})\right\}\left[1-\exp(\hbar\beta(\omega+\mathbf{k}^\intercal\mathbf{v}))\right]^{-1} \stackrel{\hbar\beta\omega\gg 1}\sim  -
\frac{2\mathrm{diag}[1,2,2]}{15 \pi\epsilon_0 }\frac{v^2}{c^2}
\frac{\omega^3}{c^3}\frac{1}{\hbar\beta\omega } $

    - $\int\frac{\mathrm{d}^2\mathbf{ k}}{(2\pi)^2}  \mathrm{Im}\left\{\underline{G}_0(\mathbf{ k}, z\to z',\omega+\mathbf{k}^\intercal\mathbf{v})\right\}\left[1-\exp(\hbar\beta(\omega+\mathbf{k}^\intercal\mathbf{v}))\right]^{-1}k_x \stackrel{\hbar\beta\omega\gg 1}\sim  -
 \frac{\mathrm{diag}\left[1,2,2\right]}{30\pi\epsilon_0}\frac{v}{c}\frac{\omega^4}{c^4}
 \frac{1}{\hbar\beta\omega} $

   * For the [plate Green's tensor](api/greenstensor?id=greenstensorplate) we find in the near-field limit ($|\mathbf{k}|^2\ll \omega_\mathrm{ch}^2/c^2 \ll \omega^2/c^2$, where $\omega_\mathrm{ch}$ is the characteristic frequency, e.g. a resonance, of the material)

    - $\underline{g}(\mathbf{k},z_a,\omega) \sim \left(\mathrm{diag}[\frac{k_x^2}{k^2},\frac{k_y^2}{k^2},1] - \underline{L}_y\right) \frac{k}{2\epsilon_0} r(\omega) e^{-2 k z_a}$ with $[\underline{L}_i]_{jk}=-\mathrm{i}\epsilon_{ijk}$ where $\epsilon_{ijk}$ is the [Levi-Civita symbol](https://en.wikipedia.org/wiki/Levi-Civita_symbol#Three_dimensions).

    - in the below  equations we will assume that $\mathrm{Im}\{r(\omega)\}\approx 2\epsilon_0\rho \omega$, with $\rho$ being a constant resistivity, and introduce $\eta = 2z_a\omega/v$ 

    - $\int\frac{\mathrm{d}^2\mathbf{ k}}{(2\pi)^2}\underline{g}_\Im(\mathbf{k},z_a,\omega+\mathbf{k}^\intercal\mathbf{v})\left[1-\exp(\hbar\beta(\omega+\mathbf{k}^\intercal\mathbf{v}))\right]^{-1} \stackrel{\hbar\beta\omega\ll 1}\sim \frac{2 v \rho}{(2z_a)^4\pi }\left(
\mathrm{diag}[\frac{\pi}{2} \eta + 4, \frac{\pi}{2} \eta + 2], \pi\eta + 6] -\underline{L}_y(\frac{3\pi}{2} + 2\eta) \right)$

    - $\int\frac{\mathrm{d}^2\mathbf{ k}}{(2\pi)^2}\underline{g}_\Im(\mathbf{k},z_a,\omega+\mathbf{k}^\intercal\mathbf{v})\frac{k_x}{1-\exp(\hbar\beta(\omega+\mathbf{k}^\intercal\mathbf{v}))} \stackrel{\hbar\beta\omega\ll 1}\sim \frac{2 v \rho}{(2z_a)^5\pi }\left(
\mathrm{diag}[\frac{9\pi}{2} + 4\eta, \frac{3\pi}{2} + 2\eta], 6\pi + 6\eta] -\underline{L}_y(\frac{3\pi}{2}\eta + 16) \right)$

    - $\int\frac{\mathrm{d}^2\mathbf{ k}}{(2\pi)^2}\underline{g}_\Im(\mathbf{k},z_a,\omega+\mathbf{k}^\intercal\mathbf{v})\left[1-\exp(\hbar\beta(\omega+\mathbf{k}^\intercal\mathbf{v}))\right]^{-1} \stackrel{\hbar\beta\omega\gg 1}\sim \frac{2 \rho}{(2z_a)^3\hbar\beta }\left(
\mathrm{diag}[1, 1, 2] -\underline{L}_y\frac{3v}{4 z_a} \right)$

## Reflection Coefficients

 - Obeys crossing relation

  * $r_\sigma(-\omega,k) = r^*_\sigma(\omega,k)$ where $\sigma=p,s$ refers to the transverse magnetic ($p$) or transverse electric ($s$) component 

 - For $k^2 \ll \omega_\mathrm{ch}^2/c^2 \ll \omega^2/c^2$ we retrieve the evanescent limit 

  * $r_p(\omega,k) \approx \frac{\epsilon(\omega) - 1}{\epsilon(\omega) +1} $ and $r_s(\omega,k) \approx (\epsilon(\omega) - 1)\frac{\omega^2}{4 c^2 k^2} + (\epsilon(\omega)^2 -1) \frac{\omega^4}{8 c^4 k^4} $

 - For $ \omega_\mathrm{ch}^2/c^2 \ll k^2 \ll \omega^2/c^2$ we retrieve the propagating limit (with respect to propagation within the material) 

  * $r_p(\omega,k) \approx \frac{\sqrt{\epsilon(\omega)} - 1}{\sqrt{\epsilon(\omega)} +1} $ and $r_s(\omega,k) \approx -\frac{\sqrt{\epsilon(\omega)} - 1}{\sqrt{\epsilon(\omega)} +1} $

## Permittivity

 - Obeys crossing relation

  * $\epsilon(-\omega) = \epsilon^*(\omega)$

## Polarizability

 - Obeys crossing relation

  * $\underline{\alpha}(-\omega) = \underline{\alpha}^*(\omega)$

 - In the case of an atom (without internal bath) moving through vacuum, we find
 
  * $\int_0^{\omega_\mathrm{cut}}\mathrm{d}\omega\,\underline{\alpha}_\Im(\omega)
\sim
\begin{cases}
\mathrm{diag}[1,1,1]\frac{\pi}{2}\alpha_0\omega_a
\quad  & \text{for}\quad \omega_\mathrm{cut}\gg \omega_a
\\\\
\alpha_0^2
  \frac{\omega_\mathrm{cut}^4}{8\pi\epsilon_0}\mathrm{diag}\left[
  \frac{ c }{3 \left(c^2-v^2\right)^2}
    ,\,
    \frac{ c  \left(c^2+v^2\right)}{3 \left(c^2-v^2\right)^3}
  ,\,
\frac{ c  \left(c^2+v^2\right)}{3 \left(c^2-v^2\right)^3}
\right]
\stackrel{v\ll c}\approx
\mathrm{diag}[1,1,1]\alpha_0^2 \frac{\omega_\mathrm{cut}^4}{24\pi\epsilon_0 c^3}
\quad  &\text{for}\quad \omega_\mathrm{cut}\ll \omega_a
\end{cases}$

 - In the case of an microscopic object with an internal bath (assuming $\mu(\omega)=\gamma=\mathrm{const.}$) we find
  
  * $\int_0^{\omega_\mathrm{cut}}\mathrm{d}\omega\,\underline{\alpha}_\Im(\omega)
\sim
\begin{cases}
\mathrm{diag}[1,1,1]
\frac{\alpha_0\omega_a^2}{2}
\left[
  \frac{\pi}{\sqrt{4 \omega_a^2-\gamma ^2}}-\frac{2 \arctan\left(\frac{\gamma ^2-2 \omega_a^2}{\gamma  \sqrt{4 \omega_a^2-\gamma ^2}}\right)}{\sqrt{4 \omega_a^2-\gamma ^2}}
\right]
\stackrel{\gamma\ll\omega_a}\sim\frac{\pi}{2}\alpha_0\omega_a
\mathrm{diag}[1,1,1]
\quad  & \text{for}\quad \omega_\mathrm{cut}\gg \omega_a
\\\\ \mathrm{diag}[1,1,1]\left(\alpha_0^2 \frac{\omega_\mathrm{cut}^4}{24\pi\epsilon_0 c^3}+ \frac{\alpha_0\gamma\omega_\mathrm{cut}^2}{2\omega_a^2}\right)\quad  &\text{for}\quad \omega_\mathrm{cut}\ll \omega_a \end{cases}$

## Memory Kernel

 - Obeys crossing relation

  * $\mu(-\omega) = \mu^*(\omega)$

## Power Spectrum

 - Power spectrum is hermitian
  * $\underline{S}^\dagger(\omega) = \underline{S}(\omega)$

 - For the static case $v\to 0$ the power spectrum takes the form
  * $\underline{S}(\omega)\sim \frac{\hbar}{\pi} \frac{\underline{\alpha}_\Im(\omega)}{1-\exp(-\hbar\beta\omega)}$

## Friction
  - For an atom (without internal bath) moving through vacuum, we can find a high temperature ($\hbar\omega_a \ll k_\mathrm{B}T$) result

   - $F\sim - \frac{v}{c} \frac{\alpha_0 }{6\pi\epsilon_0\beta} \frac{\omega_a^4}{c^4}$

  - We can retrieve asymptotical results for an atom moving above a surface for low temperatures ($\hbar\frac{v}{2z_a} \gg k_\mathrm{B}T$) as calculated [here](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.123.120401)

   - $F\sim -(63-45) \frac{\hbar \alpha_0^2\rho^2}{\pi^3} \frac{v^3}{(2z_a)^{10}} - (6-3) \frac{\alpha_0^2\rho^2}{\pi\beta^2\hbar} \frac{6v}{(2z_a)^8}$

  - Adding an Ohmic internal bath (assuming $\mu(\omega)=\gamma=\mathrm{const.}$) we find

   - $F\sim  -(63-45) \frac{\hbar \alpha_0^2\rho^2}{\pi^3} \frac{v^3}{(2z_a)^{10}} - (6-3) \frac{\alpha_0^2\rho^2}{\pi\beta^2\hbar} \frac{6v}{(2z_a)^8}-45 \frac{\hbar\alpha_0\rho}{\pi^2}\frac{\gamma}{\omega_a^2} \frac{v^3}{(2z_a)^7} - \frac{\alpha_0\rho\pi\gamma}{\hbar\beta^2\omega_a^2} \frac{8v}{(2z_a)^5}$
