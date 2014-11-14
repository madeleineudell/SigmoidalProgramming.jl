## sigmoidal functions

export Logistic, logistic, logistic_prime,
	ZeroSigmoid

type Sigmoid
	value::Function
	gradient::Function
	z::Float64
end

Logistic() = Sigmoid(logistic, logistic_prime, 0)
logistic(x) = 1/(1 + exp(-x))
logistic_prime(x) = exp(-x)/(1 + exp(-x))^2

ZeroSigmoid() = Sigmoid(x->0, x->0, 0)
