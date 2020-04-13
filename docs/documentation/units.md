# System of units
QuaCa does not anticipate any system of measurement, so essentially you have to know yourself which units you put in and which units result from this.
To make this easier, however, we **strongly** advise the use to use a system of measurement that is commonly called *natural units*, where we set $\hbar = c = k_B = 4\pi\epsilon_0 = 1$.
Please refer to the [Wikipedia page](https://en.wikipedia.org/wiki/Natural_units#%22Natural_units%22_(particle_physics_and_cosmology)) for more information.
Below you find a table which matches each parameter to it's LaTeX expressions and clarifies what we mean by it.
Further we provide a converter from common units of the respective quantity to the natural units used for QuaCa.
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
</table>
