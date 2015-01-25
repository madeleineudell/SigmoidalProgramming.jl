# SigmoidalPrograming

<!-- [![Build Status](https://travis-ci.org/madeleineudell/SigmoidalPrograming.jl.svg?branch=master)](https://travis-ci.org/madeleineudell/SigmoidalPrograming.jl) -->

SigmoidalProgramming is a Julia package for solving [sigmoidal programming problems](http://www.stanford.edu/~udell/doc/max_sum_sigmoids.pdf).
It solves problems using a branch and bound method described [here](http://www.stanford.edu/~udell/doc/max_sum_sigmoids.pdf),
and solves the convex subproblems using a cutting plane method.

Usage examples can be found in the examples directory.

## Installation

You can install the package from the Julia prompt
```
Pkg.clone("git@github.com:madeleineudell/SigmoidalProgramming.jl.git")
```