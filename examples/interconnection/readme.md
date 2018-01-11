
[paper]: https://people.eecs.berkeley.edu/~eskim/papers/HSCC18_preprint.pdf

About
============

Two examples contained in the [HSCC18 modular control systems paper][paper].

- runningmax.cc: An example to test the scalability of the interconnection decomposition. 
- consensus.cc: Multiple scalar systems are tasked with agreeing on a consensus location. The only global information available to each individual system is the average state amongst all systems. 


Installation Instructions 
============

1. To run these examples, you need the [CUDD](http://vlsi.colorado.edu/~fabio/CUDD/) library installed. On MacOS and using [Homebrew](https://brew.sh/) you can run 
```
brew install cudd 
```
in a terminal then
```
brew info cudd
```
to find the library location.

2. Edit the Makefile 
   
   * Set CC to a compiler that supports C++11, e.g. `CC = clang++`
   
   * set CUDDPATH to the location of the CUDD library, e.g. `CUDDPATH = /usr/local/Cellar/cudd/3.0.0`

3. Compile `runningmax.cc` and `consensus.cc` by running

  ```
  $ make 
  ``` 

4. Execute examples with

  ```
  $./runningmax.o
  $./consensus.o
  ```


5. Visualize consensus example
  
  ```
  python visualize_consensus.py
  ```