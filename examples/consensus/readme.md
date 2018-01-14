
[paper]: https://people.eecs.berkeley.edu/~eskim/papers/HSCC18_preprint.pdf

About
============

Consensus example contained in the [HSCC18 modular control systems paper][paper].

* consensus.cc: Multiple scalar systems are tasked with agreeing on a consensus location. The only global information available to each individual system is the average state amongst all systems. 


Installation Instructions 
============

1. To run these examples, you need the [CUDD](http://vlsi.colorado.edu/~fabio/CUDD/) library installed. On MacOS and using [Homebrew](https://brew.sh/) you can run `brew install cudd`
in a terminal then `brew info cudd` to find the library location.

2. Edit the Makefile 
  
  * Set CC to a compiler that supports C++11, e.g. `CC = clang++`
  
  * set CUDDPATH to the location of the CUDD library, e.g. `CUDDPATH = /usr/local/Cellar/cudd/3.0.0`

How to Run
=============

1. Compile `consensus.cc` by running `$ make`

2. Execute examples with `$./consensus.o`

3. To recreate Figure 7, run `python visualize_consensus.py`. This requires the python graphical plotting library matplotlib.