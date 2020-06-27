## sigmoidal functions

export logistic, logistic_prime, weibull, weibull_prime

logistic(x) = 1/(1 + exp(-x))
logistic_prime(x) = exp(-x)/(1 + exp(-x))^2
weibull(x, coef, scale, shape) = (1 - exp(-((x/scale)^shape)))*coef
weibull_prime(x, coef, scale, shape) = coef*(shape/scale)*((x/scale)^(shape-1))*exp(-((x/scale)^shape))
