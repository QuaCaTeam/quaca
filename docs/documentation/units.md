# System of units
QuaCa does not anticipate any system of measurement, so essentially you have to know yourself which units you put in and which units result from this.
To make this easier, however, we **strongly** advise the use to use a system of measurement that is commonly called *natural units*, where we set $\hbar = c = k_B = 1$.
Please refer to the [Wikipedia page](https://en.wikipedia.org/wiki/Natural_units#%22Natural_units%22_(particle_physics_and_cosmology)) for more information.
With this small tool you can convert from more common SI units to the natural units needed for QuaCa

<table>
  <tr>
    <th>Input value</th>
    <th>Input unit</th>
    <th>Value for QuaCa</th>
  </tr>
  <tr>
    <td><input type="text" placeholder="value" oninput="unitConverter()" onchange="unitConverter()" size="10" id="inputValue"></td>
    <td>
    <select id="inputUnit" onchange="unitConverter()">
    <option value="eV">eV</option>
    <option value="nm">nm</option>
    </select>
    </td>
    <td><p><span id="outputValue">0</span></p></td>
  </tr>
</table>


# Reference table
QuaCa tries to write its variable names just like in LaTeX.
Below you'll find a table which matches each parameter to it's LaTeX expressions and clarifies what we mean by it.


| Code                | Math       | Physical quantity             |
|---------------------|------------|-------------------------------|
| `double omega_a`    | $\omega_a$ | Resonance frequency of dipole |
| `double alpha_zero` | $\alpha_0$ | Vacuum polarizability         |
| `double gamma`      |            |                               |
