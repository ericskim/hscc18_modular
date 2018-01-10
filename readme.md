# Modular SCOTS

[paper]: https://people.eecs.berkeley.edu/~eskim/papers/HSCC18_preprint.pdf
> This repository is a submission to the [HSCC2018 Repeatability Evaluation](https://www.hscc2018.deib.polimi.it/repeatability-evaluation).

This is a modified version of [SCOTSv0.2](https://gitlab.lrz.de/matthias/SCOTSv0.2), which is an open source software tool to compute discrete abstractions and symbolic controllers. Whereas the vanilla implementation treats control systems as monolithic objects, this version places an emphasis on computing discrete abstractions of individual subsystems, then interconnecting them together. Examples highlight how this reduces the controller abstraction time. The underlying theory behind this approach can be found in [this paper][paper].

Bug reports and feature requests can be submitted to <eskim@eecs.berkeley.edu> 

# Additional Functionality 

Most of the functionality of the original version has been retained, but some features to take advantage of the modular approach to constructing control system abstractions can be found in the files 

- [./src/FunctionAbstracter.hh](./src/FunctionAbstracter.hh)
- [./src/EnfPre.hh](./src/EnfPre.hh)

Two examples from the paper and instructions to run them can be found in the  [./examples/interconnection/](./examples/interconnection/) directory.

For implementation details please have a look in the C++ documentation ./doc/html/index.html

### How to use:

* The basic implementation of **SCOTS** is inlined and header only. Hence, only a working C++ compiler
  with C++11 support is needed.

* The best way to find out about **SCOTS** is to clone the repository 
  
    `$ git clone https://gitlab.lrz.de/matthias/SCOTSv0.2.git`
  
    and run some of the examples in the example directory: 

  * ./examples/dcdc/
  * ./examples/vehicle/

    Have a look in the readme file for some info and compiler options
