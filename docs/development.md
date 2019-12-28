# Project structure

## Folder structure

The project is organized as follows:

```bash
    QuaCa/
    ├── app
    ├── bin
    ├── build
    ├── data
    ├── docs
    ├── include
    ├── lib
    ├── src
    └── test
```

The directory contains ...

- `app/`: source files for the executables
- `bin/`: the executables (binary files)
- `build/`: build files produced by cmake
- `data/`: various data files for the tests, tutorials and other applications
- `docs/`: the documentation
- `include/`: library headers
- `src/`: source files for the quaca library
- `test/`: unit tests

# Current tasks

```mermaid
gantt
       dateFormat  DD-MM-YY
       title Development tasks

       section Marty
       Testliste            :done,    des1, 01-12-19, 17-12-19
       Asymptoten ausrechnen :active,  des1, 17-12-19, 3w
       
       section Simon
       Integrationsroutine            :done,    des1, 01-12-19, 17-12-19
       Tests fuer Greenstensor        :active,  des1, 17-12-19, 3w
       
       section Christoph
       Dokumentation aufsetzen  :done,    des1, 01-12-19, 17-12-19
       Polarisierbarkeit        :active,  des1, 17-12-19, 3w
```
