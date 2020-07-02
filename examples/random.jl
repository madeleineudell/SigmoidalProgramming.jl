## Let's solve a random problem
using SigmoidalProgramming
using Random
import DataStructures: dequeue!

# problem data
Random.seed!(4)
nvar = 200
nconstr = 40

l = -rand(nvar)
u = rand(nvar)
A = rand(nconstr, nvar)
b = rand(nconstr)
z = zeros(nvar)
fs = fill(logistic, nvar)
dfs = fill(logistic_prime, nvar)
problem = LinearSP(fs, dfs, z, A, b)

# branch and bound
pq, bestnodes, lbs, ubs = solve_sp(l, u, problem)

node = dequeue!(pq)
# println("best node has node.ub = $(node.ub) and solution $(node.x)")
println("lbs: ",lbs)
println("ubs: ",ubs)
