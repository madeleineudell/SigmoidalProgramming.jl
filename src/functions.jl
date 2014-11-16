## sigmoidal functions
using Convex

export Logistic, logistic, logistic_prime,
	ZeroSigmoid, SigmoidalFunction

type SigmoidalFunction
	value::Function
	gradient::Function
	z::Float64
end

Logistic() = SigmoidalFunction(logistic, logistic_prime, 0)
logistic(x) = 1/(1 + exp(-x))
logistic_prime(x) = exp(-x)/(1 + exp(-x))^2

ZeroSigmoid() = SigmoidalFunction(x->0, x->0, 0)

function SigmoidalFunction(x::AbstractExpr)
	size(x) == (1,1) && length(x.children) == 1 && size(x.children[1]) == (1,1) ||
		error("SigmoidalProgramming requires every term added to a sum-of-sigmoids 
			   objective to be univariate, with univariate argument")
	copy_of_x = typeof(x)(Variable(x.children[1].size))
	if vexity(x) == ConvexVexity() || AffineVexity()
		f = SigmoidalFunction(v->(copy_of_x.children[1].value = v; evaluate(copy_of_x)),
			                  v->0, # slope of a convex function should never be used
			                  float("Inf"))
	elseif vexity(x) == ConcaveVexity()
		error("SigmoidalProgramming does not currently support concave objective terms")
	end
	return f
end