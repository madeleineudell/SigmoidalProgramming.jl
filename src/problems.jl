export SigmoidalProgram, LinearSP

## problem types
abstract type SigmoidalProgram end

# LinearSP represents the problem
# maximize    sum_i f_i(x_i)
# subject to  A x <= b
#             C x == d
struct LinearSP <: SigmoidalProgram
    fs::Array{Function,1}
    dfs::Array{Function,1}
    z::Array{Float64,1}
    A::Array{Float64,2}
    b::Array{Float64,1}
    C::Array{Float64,2}
    d::Array{Float64,1}
    function LinearSP(fs,dfs,z,A,b,C,d)
        nconstr, nvar = size(A)
        @assert(length(b) == nconstr)
        nconstr, nvar = size(C)
        @assert(length(d) == nconstr)
        @assert(length(fs) == nvar)
        @assert(length(dfs) == nvar)
        @assert(length(z) == nvar)
        new(fs,dfs,z,A,b,C,d)
    end
end
# Default constructor when no equality constraints are present
LinearSP(fs, dfs, #::Array{Function,1},dfs::Array{Function,1},
         z::Array{Float64,1},A::Array{Float64,2},b::Array{Float64,1}) =
         LinearSP(fs,dfs,z,A,b,zeros(0,length(fs)),zeros(0))

function addConstraints!(m::Model, x, p::LinearSP)
    # Ax <= b
    nconstr, nvar = size(p.A)
    for i=1:nconstr
        @constraint(m, sum(p.A[i,j]*x[j] for j=1:nvar) <= p.b[i])
    end
    # Cx == d
    nconstr, nvar = size(p.C)
    for i=1:nconstr
        @constraint(m, sum(p.C[i,j]*x[j] for j=1:nvar) == p.d[i])
    end
end
