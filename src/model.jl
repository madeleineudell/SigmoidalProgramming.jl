using Convex

import Base.sign
export SigmoidalVexity, SumofSigmoidalVexity,
	SigmoidalAtom, evaluate, sign, monotonicity, curvature, conic_form!,
	logisticsigmoid, zerosigmoid,
	solve_sp!

type SigmoidalVexity <: Vexity            end
type SumofSigmoidalVexity <: Vexity            end

+(s::SigmoidalVexity, t::SigmoidalVexity) = SumofSigmoidalVexity()
+(s::SigmoidalVexity, t::ConvexVexity) = SumofSigmoidalVexity()
+(s::SigmoidalVexity, t::ConcaveVexity) = SumofSigmoidalVexity()
+(s::SigmoidalVexity, t::SumofSigmoidalVexity) = SumofSigmoidalVexity()

type SigmoidalAtom <: AbstractExpr
  head::Symbol
  sigmoid::SigmoidalFunction
  id_hash::Uint64
  children::(AbstractExpr,)
  size::(Int64, Int64)

  function SigmoidalAtom(s::SigmoidalFunction, x::AbstractExpr)
    if x.size!=(1,1)
    	error("The argument to a sigmoid must be univariate; got $x")
    end
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

logisticsigmoid(x::AbstractExpr) = SigmoidalAtom(Logistic(), x)
zerosigmoid(x::AbstractExpr) = SigmoidalAtom(ZeroSigmoid(), x)

## Solve

function solve_sp!(problem::Problem, l, u; kwargs...)

	# check we have a sigmoidal program
	if !(vexity(problem.objective) == SumofSigmoidalVexity()) || 
		!all(Bool[vexity(c) == AffineVexity() for c in problem.constraints])
			error("Not a sigmoidal program")
	end
	if problem.sense == :minimize
		error("only sigmoidal maximization is supported for now")
	end
	
	# children gives the expressions corresponding to the arguments to the sigmoids
	sigmoids, children, constraints = gather_sigmoids(problem.objective)
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

function gather_sigmoids(obj::AbstractExpr)
	if obj.head == :sigmoid
		return [obj.sigmoid], [obj.children...], Constraint[]
	elseif vexity(obj) == Convex() || vexity(obj) == Affine()
		f = SigmoidalFunction(obj)
		return [f], [obj.children...], Constraint[]
	elseif vexity(obj) == Concave()
		# add a new hypograph variable
	elseif obj.head == :sum || obj.head == :+
		s = SigmoidalFunction[]
		c = AbstractExpr[]
		constr = Constraint[]
		for child in obj.children
			si,ci,constri = gather_sigmoids(child) 
			append!(s,si)
			append!(c,ci)
			append!(constr,constri)
		end
		return s,c,constr
	else
		error("$obj is not a sum of sigmoids")
	end
end