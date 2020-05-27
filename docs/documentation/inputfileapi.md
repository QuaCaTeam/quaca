# Input File {docsify-ignore-all}

As already indicated in [Documentation](api/greenstensor), we have to provide proper input files for the calculations of QuaCa. In the following, we present some useful input files which yield already published data.

## Atom above a Drude Bulk
This configuration yields the numerical results given in [this publication](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.123.120401)
``` json
{
    "GreensTensor": {
        "type": "plate",
        "v": 1e-4,
        "beta": 1e6,
        "za": 0.01,
        "delta_cut": 20,
        "rel_err_0": 1e-4,
        "rel_err_1": 1e-2
    },
    "ReflectionCoefficients": {
        "type": "local bulk"
    },
    "Permittivity": {
        "type": "drude",
        "omega_p": 9.0,
        "gamma": 0.1
    },
    "Polarizability": {
        "omega_a": 1.3,
        "alpha_zero": 6e-9
    },
    "Friction": {
        "relerr_omega": 1e-1
    },
    "Looper": {
        "type": "v",
        "scale": "log",
        "start": 1e-4,
        "end": 1e-2,
        "steps": 40
    }
}
```
