# Modular SCOTS

[paper]: https://people.eecs.berkeley.edu/~eskim/papers/HSCC18_preprint.pdf
> This repository is a submission to the [HSCC2018 Repeatability Evaluation](https://www.hscc2018.deib.polimi.it/repeatability-evaluation).

This is a modified version of [SCOTSv0.2](https://gitlab.lrz.de/matthias/SCOTSv0.2), which is an open source software tool for controller synthesis by first abstracting systems into discrete automata. Whereas the regular implementation treats control systems as monolithic objects, this version places an emphasis on computing discrete abstractions of individual subsystems, then interconnecting them together. Examples highlight how this approach reduces the control system abstraction time and the underlying theory can be found in [this paper][paper].

Bug reports and feature requests can be submitted to <eskim@eecs.berkeley.edu> 

# Additional Functionality 

Most of the functionality of the original version has been retained, but features to take advantage of the modular approach to constructing control system abstractions can be found in the files 

- [./src/FunctionAbstracter.hh](./src/FunctionAbstracter.hh) for abstracting both control systems and interconnection maps. 
- [./src/EnfPre.hh](./src/EnfPre.hh) for using computing controllable predecessors. Predecessors are used in controller synthesis via fixed point computations. 

For implementation details please have a look in the C++ [documentation](./doc/html/index.html). 

### How to use:

* The basic implementation of **SCOTS** is inlined and header only. Hence, only a working C++ compiler
  with C++11 support is needed.

* The best way to find out about **SCOTS** is to clone the repository and run the examples

  - [./examples/dcdc_bdd/](./examples/dcdc_bdd/)
  - [./examples/vehicle_bdd/](./examples/vehicle_bdd)

* Two additional examples from the [paper][paper] and instructions on how to run them can be found in the  [./examples/interconnection/](./examples/interconnection/) directory.


