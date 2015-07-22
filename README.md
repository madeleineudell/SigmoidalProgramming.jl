# SigmoidalPrograming

<!-- [![Build Status](https://travis-ci.org/madeleineudell/SigmoidalPrograming.jl.svg?branch=master)](https://travis-ci.org/madeleineudell/SigmoidalPrograming.jl) -->

SigmoidalProgramming is a Julia package for solving [sigmoidal programming problems](http://www.stanford.edu/~udell/doc/max_sum_sigmoids.pdf).
It solves problems using a branch and bound method described [here](http://www.stanford.edu/~udell/doc/max_sum_sigmoids.pdf),
and solves the convex subproblems using a cutting plane method.

Usage examples can be found in the examples directory.

## Installation

You can install the package from the Julia prompt
```
Pkg.add("SigmoidalProgramming")
```

## Linear Sigmoidal Programs

The `LinearSP` problem class can be used to solve any sigmoidal programming problem
with linear equality and inequality constraints, of the form
```
    minimize       sum_i f_i(x_i)
    subject to     A x <= b
                   C x == d
                   l <= x <= u
```
with variable x, where `f_i` is a sigmoidal function for every `i`. 
A *sigmoidal* function is a function `f(x)` that, for some number `z`, is convex 
in `x` for `x <= z` and is concave in `x` for `x >= z`.
(Convex functions are sigmoidal with `z = float("inf")`, and concave functions are sigmoidal with `z = -float("inf")`.)

The functions `f` are specified with three lists:

1. `fs`, a list of functions such that `fs[i](x)` computes the value of `f_i` at `x`.  
2. `dfs`, a list of functions such that `dfs[i](x)` computes the derivative of `f_i` at `x`.  
3. `z`, a list of numbers such that `f_i` is convex 
in its argument `x` for `x <= z` and is concave in `x` for `x >= z`.

Let's see how to use it to solve a random sigmoidal program with only inequality constraints.
We'll let `f_i = x -> 1/(1 + exp(-x))` be the `logistic` function for every `i`. 
The logistic function and its derivative `logistic_prime` are defined in the file `functions.jl`.
```
using SigmoidalProgramming

# generate problem data
srand(4)
nvar = 200
nineqconstr = 20
l = -rand(nvar)
u = rand(nvar)
A = rand(nineqconstr, nvar)
b = rand(nineqconstr)
z = zeros(nvar)
fs = fill(logistic, nvar)
dfs = fill(logistic_prime, nvar)
problem = LinearSP(fs, dfs, z, A, b)

# branch and bound to solve the problem
# pq is a priority queue of the branch and bound nodes at the leaves of the tree
# bestnodes is a list of the best branch and bound nodes found, in the order they were found
pq, bestnodes, lbs, ubs = solve_sp(l, u, problem)

# the best node found yet is the top node on the priority queue
node = dequeue!(pq)
# println("best node has node.ub = $(node.ub) and solution $(node.x)")

# lbs and ubs record the upper and lower bounds on the optimal value
# found at each iteration
println("lbs: ",lbs)
println("ubs: ",ubs)
```

A problem with equality constraints `C x = d` as well can be solved using `problem = LinearSP(fs, dfs, z, A, b, C, d)
`

A few simple examples of usage can be found in the examples directory.

## Solver parameters

The accuracy of the solution can be controlled with the parameters

* `maxiters`, the maximum number of branch and bound iterations (default: `100`)
* `TOL`, the desired gap between the upper and lower bound on the optimal value (default: `1e-2`)

The solver terminates as soon as either `maxiters` branch and bound nodes have been split,
or the gap between the upper and lower bound on the optimal value has been proved to be less
than `TOL`.

The call signature for `solve_sp` is as follows:
```
solve_sp(l, u, problem::SigmoidalProgram; 
                  TOL = 1e-2, maxiters = 100, verbose = false)
```

