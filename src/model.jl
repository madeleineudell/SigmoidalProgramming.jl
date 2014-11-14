using Convex

type SigmoidalVexity <: Vexity            end
type SumofSigmoidalVexity <: Vexity            end

+(s::SigmoidalVexity, t::SigmoidalVexity) = SumofSigmoidalVexity()
+(s::SigmoidalVexity, t::SumofSigmoidalVexity) = SumofSigmoidalVexity()

type SigmoidalAtom <: AbstractExpr
  head::Symbol
  sigmoid::Sigmoid
  id_hash::Uint64
  children::(AbstractExpr,)
  size::(Int64, Int64)

  function SigmoidalAtom(s::Sigmoid, x::AbstractExpr)
    if x.size!=(1,1)
    	error("The argument to a sigmoid must be univariate; got $x")
    children = (x,)
    return new(:sigmoid, s, hash(children), children, (1,1))
  end
end

function sign(x::SigmoidalAtom)
  return NoSign()
end

function monotonicity(x::SigmoidalAtom)
  return NotMonotonic()
end

function curvature(x::SigmoidalAtom)
  return NotDcp()
end

function evaluate(x::SigmoidalAtom)
  return x.sigmoid.value(x.children[1])
end

function conic_form!(x::SigmoidalAtom, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, x)
    # we just need to add the constraints for the children to unique_conic_forms
    c = x.children[1]
    conic_form!(c, unique_conic_forms)
  end
  return get_conic_form(unique_conic_forms, c)
end

## Particular sigmoids

logistic(x::AbstractExpr) = SigmoidalAtom(Logistic(), x)
zerosigmoid(x::AbstractExpr) = SigmoidalAtom(ZeroSigmoid(), x)

## Solve

function solve_sp(problem::Problem, l, u; kwargs...)

	# check we have a sigmoidal program
	if vexity(problem.objective) == SumofSigmoidalVexity() && 
		all(Bool[vexity(c) == AffineVexity() for c in problem.constraints])
		continue
	else
		error("Not a sigmoidal program")
	end

	# children gives the expressions corresponding to the arguments to the sigmoids
	sigmoids, children = gather_sigmoids(problem.objective)
	problem.objective = Constant(0)

	c, A, b, cones, var_to_ranges, vartypes = conic_problem(problem)

	full_sigmoids = fill(ZeroSigmoid(), length(c))
	for i=1:length(sigmoids)
		child_idx = var_to_ranges[children[i]][1]
		full_sigmoids[child_idx] = sigmoids[child_idx]
	end

	# unfold linear equality and inequality constraints
	sp = LinearSP(full_sigmoids,A,b,cones)

	pq, bestnodes, lbs, ubs = 
	              solve_sp(l, u, sp; kwargs...)

	problem.solution = Solution(bestnodes[end].x, :Optimal, lbs[end])
	populate_variables!(problem, var_to_ranges)

	# Populate the problem with the solution
	problem.optval = lbs[end]
	problem.status = problem.solution.status
end
