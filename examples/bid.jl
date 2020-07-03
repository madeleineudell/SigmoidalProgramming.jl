## Let's solve the bidding problem
using SigmoidalProgramming
using Random

# problem data
Random.seed!(4)

function bidding(nvar=500; problemtype=:paper)
    nconstr = 1
    l = fill(-10, nvar)
    u = fill(10, nvar)
    A = ones(nconstr, nvar)

    # boring bidding
    if problemtype==:boring
    B = [.4*nvar]
    z = zeros(nvar)
    fs = Function[x -> i/(1 + exp(-x)) for i=1:nvar]
    dfs = Function[x-> i*exp(-x)/(1 + exp(-x))^2 for i=1:nvar]

    # bidding example from the paper
    elseif problemtype==:paper
        v = 4*rand(nvar)
        beta = -3*v
        alpha = 10
        B = [.2*sum(v)]
        z = -beta/alpha
        fs = Function[x->logistic(alpha*x + beta[i]) - logistic(beta[i]) for i=1:nvar]
        dfs = Function[x->alpha*logistic_prime(alpha*x + beta[i]) for i=1:nvar]
    else
        error("did not recognize problem type $problemtype")
    end

    problem = LinearSP(fs, dfs, z, A, B)

    # branch and bound
    pq, bestnodes, lbs, ubs = @time solve_sp(l, u, problem, verbose=false;
                                             TOL=nvar*1e-3, maxiters=1000)
    println("$(length(lbs)) iterations")
    node = bestnodes[end] #dequeue!(pq)
    # println("best node has node.ub = $(node.ub) and solution $(node.x)")
    println("lb: ",lbs[end])
    println("ub: ",ubs[end])
    println("ratio: ",(ubs[end]-lbs[end])/lbs[end])

    return pq, bestnodes, lbs, ubs
end

# for size=[1,10,20,50,100,500,1000,2000,5000,10000]
#     println("bidding $size")
#     bidding(size)
# end

pq, bestnodes, lbs, ubs = bidding(20)
