# Reference table
QuaCa tries to write its variable names just like in LaTeX.
Below you'll find a table which matches each parameter to it's LaTeX expressions and clarifies what we mean by it.


| Code                | Math       | Physical quantity             |
|---------------------|------------|-------------------------------|
| `double omega_a`    | $\omega_a$ | Resonance frequency of dipole |
| `double alpha_zero` | $\alpha_0$ | Vacuum polarizability         |
| `double gamma`      |            |                               |


# Conversion table


<table>
  <tr>
    <th>Math</th>
    <th></th>
    <th>Code</th>
  </tr>
  <tr>
    <td>$\omega_a$ = <input type="text" placeholder="frequency" oninput="frequencyConverter(this.value)" onchange="frequencyConverter(this.value)" size="10"></td>
    <td>
    <select>
    <option value="eV">eV</option>
    <option value="nm">nm</option>
    </select>
    </td>
    <td><p><code>omega_a</code> = <span id="outputFrequency">0</span></p></td>
  </tr>
</table>
