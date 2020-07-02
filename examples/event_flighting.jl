using SigmoidalProgramming
using Dates, CSV, Random, DataFrames

export event_flighting


"to find the optimal flighting for a single event"
function event_flighting(
    scale::Float64,
    shape::Float64,
    coef::Float64,
    retention::Float64,
    spend::Float64;
    verbose::Int=0,
    maxiters::Int64,
    TOL::Float64=0.01,
    nweeks::Int=52,
    nconstr::Int=1,
    lower_budget::Float64=0.5,
    upper_budget::Float64=1.5,
    n_segments::Int=20,
)

    l = fill(0, 3 * nweeks)
    u = cat(fill(spend, 2 * nweeks), fill(spend / 5, nweeks), dims=(1,))
    A = hcat(zeros(nconstr, nweeks * 2), ones(nconstr, nweeks))

    z = cat(fill(scale*((shape - 1) / shape) ^ (1 / shape), nweeks * 2), fill(0, nweeks), dims=(1,))

    zz = zeros(nweeks * 3)
    C = copy(zz)
    C[1] = 1
    C[105] = -1
    C = transpose(C)
    D = fill(0, nweeks * 2)

    for i = 2 : nweeks
        c = copy(zz)
        c[i-1] = -retention
        c[i] = 1
        c[i + nweeks * 2] = -1
        C = vcat(C, transpose(c))
    end

    for i = nweeks + 1 : nweeks * 2
        c = copy(zz)
        c[i - 1] = -retention
        c[i] = 1
        C = vcat(C, transpose(c))
    end

    fs1 = Function[x -> weibull(x, coef, scale, shape) for i=1 : nweeks * 2]
    fs2 = Function[x -> 0 for i=1 : nweeks]
    fs = vcat(fs1, fs2)
    dfs1 = Function[x -> weibull_prime(x, coef, scale, shape) for i=1 : nweeks * 2]
    dfs2 = Function[x -> 0 for i=1 : nweeks]
    dfs = vcat(dfs1, dfs2)

    output_curve = DataFrame()
    iterlog = DataFrame()

    for B = range(
        spend * lower_budget,
        spend * upper_budget,
        step=spend * (upper_budget - lower_budget) / n_segments
    )
        println("budget: ", B)
        println(now())

        problem = LinearSP(fs, dfs, z, A, [B], C, D)

        l = fill(0, 3 * nweeks)
        u = cat(fill(spend, 2 * nweeks), fill(spend / 2, nweeks), dims=(1,))

        # branch and bound
        pq, bestnodes, lbs, ubs, log = @time solve_sp(
            l, u, problem; TOL=TOL, maxiters=maxiters, verbose=verbose
        )

        grps = DataFrame(
            period = 1 : nweeks,
            spend=fill(B, nweeks),
            grps=bestnodes[end].x[nweeks * 2 + 1: nweeks * 3],
            lb=fill(lbs[end], nweeks)
        )
        output_curve = vcat(output_curve, grps)
        log[!, :budget] .= B
        iterlog = vcat(iterlog, log)
    end

    println(now())
    return iterlog, output_curve
end
