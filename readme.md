# Welcome to SCOTSv0.2!

> SCOTSv0.2 is currently under development. Although the core
> algorithms are stable and usable, the user interfaces might change...


**SCOTS** is an open source software tool to compute discrete abstractions and symbolic controllers!

Please consult the [manual](https://gitlab.lrz.de/matthias/SCOTSv0.2/raw/master/manual/manual.pdf) for installation instructions,
usage description and background information.

For implementation details please have a look in the C++ documentation ./doc/html/index.html

Bug reports and feature requests are happily received at <matthias.rungger@tum.de> 

### How to use:

* The basic implementation of **SCOTS** is inlined and header only. Hence, only a working C++ compiler
  with C++11 support is needed.

* The best way to find out about **SCOTS** is to clone the repository 
  
    `$ git clone https://gitlab.lrz.de/matthias/SCOTSv0.2.git`
  
    and run some of the examples in the example directory: 

  * ./examples/dcdc/
  * ./examples/vehicle/

    Have a look in the readme file for some info and compiler options
  
### What's new:

* New data structure to store the transition function of symbolic models
   
* New synthesis algorithms for invariance and reachability specifications 
    (see the manual for details)

* Dynamic variable reordering can now be safely activated throughout all computations in the BDD implementation

* An example demonstrating the usage of validated ODE solvers 
    to compute a priori enclosures and growth bounds

* Complete redesign of the code base to accommodate for modern C++

* Doxygen documentation in ./doc/html/

