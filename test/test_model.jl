using SigmoidalProgramming
using Convex

# problem data
nvar = 5
nconstr = 1
l = zeros(nvar)
u = ones(nvar)
A = ones(nconstr, nvar)
b = nvar/2*ones(nconstr)

x = Variable(nvar)
p = maximize(sum([logisticsigmoid(x[i]) for i=1:nvar]), A*x <= b)
solve_sp!(p, l, u)

@show p.optval, x.value