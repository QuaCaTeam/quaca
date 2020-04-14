# Input {docsify-ignore-all}

In order to lift the curse of looking up the conversion tables to correctly implement your data into QuaCa, in the following we provide a short note on the used system of units, including a useful converter, as well as mention which quantities are not burdened by units.

# System of Units
In QuaCa the used system of units are the so-called *natural units*, where we set $\hbar = c = k_B = 4\pi\epsilon_0 = 1$.
Please refer to the [Wikipedia page](https://en.wikipedia.org/wiki/Natural_units#%22Natural_units%22_(particle_physics_and_cosmology)) for more information.
Below you find a table which matches each parameter to it's LaTeX expressions and clarifies its meaning.
Further we provide a converter from several common units of the respective quantity of the natural units used for QuaCa.
Quantities without converter in its row, have the same unit as the quantity above.

<table>
  <tr>
    <th>Code</th>
    <th>Math</th>
    <th>Physical quantity</th>
    <th>Input value</th>
    <th>Input unit</th>
    <th>Value for QuaCa</th>
  </tr>
  <tr>
    <td><code>double v</code></td>
    <td>$v$</td>
    <td>Velocity of the atom</td>
    <td><input type="text" placeholder="value" oninput="unitConverterVel()" onchange="unitConverterVel()" size="10" id="inputValueVel"></input></td>
    <td>
    <select id="inputUnitVel" onchange="unitConverterVel()">
    <option value="m/s">m/s</option>
    <option value="c">c</option>
    </select>
    </td>
    <td><p><span id="outputValueVel">0</span></p></td>
  </tr>
  <tr>
    <td><code>double beta</code></td>
    <td>$\beta$</td>
    <td>Inverse temperature</td>
    <td><input type="text" placeholder="value" oninput="unitConverterTemp()" onchange="unitConverterTemp()" size="10" id="inputValueTemp"></input></td>
    <td>
    <select id="inputUnitTemp" onchange="unitConverterTemp()">
    <option value="1/K">1/K</option>
    <option value="eV-1">eV^-1</option>
    </select>
    </td>
    <td><p><span id="outputValueTemp">0</span></p></td>
  </tr>
  <tr>
    <td><code>double omega_a</code></td>
    <td>$\omega_a$</td>
    <td>Resonance frequency of atom</td>
    <td><input type="text" placeholder="value" oninput="unitConverterFreq()" onchange="unitConverterFreq()" size="10" id="inputValueFreq"></input></td>
    <td>
    <select id="inputUnitFreq" onchange="unitConverterFreq()">
    <option value="Hz">Hz</option>
    <option value="eV">eV</option>
    </select>
    </td>
    <td><p><span id="outputValueFreq">0</span></p></td>
  </tr>
  <tr>
    <td><code>double omega_0</code></td>
    <td>$\omega_0$</td>
    <td>Central frequency of the material</td>
    <td></td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td><code>double omega_p</code></td>
    <td>$\omega_p$</td>
    <td>Plasma frequency of the material</td>
    <td></td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td><code>double gamma</code></td>
    <td>$\gamma$</td>
    <td>Damping constant</td>
    <td></td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td><code>double alpha_zero</code></td>
    <td>$\alpha_0$</td>
    <td>Static polarizability of the atom</td>
    <td><input type="text" placeholder="value" oninput="unitConverterPol()" onchange="unitConverterPol()" size="10" id="inputValuePol"></input></td>
    <td>
    <select id="inputUnitPol" onchange="unitConverterPol()">
    <option value="ang3">4pi*eps0*Ang^3</option>
    <option value="au">a.u.</option>
    <option value="eV-3">eV^-3</option>
    </select>
    </td>
    <td><p><span id="outputValuePol">0</span></p></td>
  </tr>
  <tr>
    <td><code>double za</code></td>
    <td>$z_a$</td>
    <td>Distance between atom and surface</td>
    <td><input type="text" placeholder="value" oninput="unitConverterLen()" onchange="unitConverterLen()" size="10" id="inputValueLen"></input></td>
    <td>
    <select id="inputUnitLen" onchange="unitConverterLen()">
    <option value="m">m</option>
    <option value="eV-1">eV^-1</option>
    </select>
    </td>
    <td><p><span id="outputValueLen">0</span></p></td>
  </tr>
  <tr>
    <td><code>double thickness</code></td>
    <td>$d$</td>
    <td>Finite thickness of the material slab</td>
    <td></td>
    <td></td>
    <td></td>
  </tr>
</table>

## Unitless Quantities

Furthermore, the code contains the following unitless quantities
<table>
  <tr>
    <th>Code</th>
    <th>Math</th>
    <th>Physical quantity</th>
  </tr>
  <tr>
    <td><code>double eps_inf</code></td>
    <td>$\epsilon_\infty$</td>
    <td>High-frequency permittivity limit</td>
  </tr>
  <tr>
    <td><code>double delta_cut</code></td>
    <td>$\delta_\mathrm{cut}$</td>
    <td>Numerical cut-off of the real part $\kappa$ integration</td>
  </tr>
  <tr>
    <td><code>double rel_err_0</code></td>
    <td>$\delta_\mathrm{rel}^0$</td>
    <td>Relative accuracy of the first (deepest) integration of the Green's Tensor</td>
  </tr>
  <tr>
    <td><code>double rel_err_1</code></td>
    <td>$\delta_\mathrm{rel}^1$</td>
    <td>Relative accuracy of the second (deepest) integration of the Green's Tensor</td>
  </tr>
  <tr>
    <td><code>double relerr_omega</code></td>
    <td>$\delta_\omega$</td>
    <td>Relative accuracy of the $\omega$ integration</td>
  </tr>
</table>
