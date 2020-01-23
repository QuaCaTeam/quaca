# System of units
We insert the parameters needed for the calculation of quantum friction in an system of units that is commonly called *natural units*.
So we set $\hbar = c = 1$.

# Reference table
QuaCa tries to write its variable names just like in LaTeX.
Below you'll find a table which matches each parameter to it's LaTeX expressions and clarifies what we mean by it.


| Code                | Math       | Physical quantity             |
|---------------------|------------|-------------------------------|
| `double omega_a`    | $\omega_a$ | Resonance frequency of dipole |
| `double alpha_zero` | $\alpha_0$ | Vacuum polarizability         |
| `double gamma`      |            |                               |


# Conversion table
Below you can convert the parameters from more common units into ones QuaCa demands.

<table>
  <tr>
    <th>Math</th>
    <th></th>
    <th>Code</th>
  </tr>
  <tr>
    <td>$\omega_a$ = <input type="text" placeholder="frequency" oninput="frequencyConverter()" onchange="frequencyConverter()" size="10" id="inputFrequency"></td>
    <td>
    <select id="inputFrequencyUnit" onchange="frequencyConverter()">
    <option value="eV">eV</option>
    <option value="nm">nm</option>
    </select>
    </td>
    <td><p><code>omega_a</code> = <span id="outputFrequency">0</span></p></td>
  </tr>
</table>
