# Input File {docsify-ignore-all}

The provided ready to use application to compute either the friction force or the decay rate require to provide a json formated input file. Furthermore, all the individual classes can be instantiated with a given input field.

The general structure of the input file is such, that any class instantiated within any of the application, gets a section in the input file, where all the relevant parameters are defined. The avaiable parameters for each class can be seen in the [API documentation](../api/greenstensor.md) under the section input file.
In the following, we present an example input file, which yields already published data.

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
