## sigmoidal functions

export logistic, logistic_prime

logistic(x) = 1/(1 + exp(-x))
logistic_prime(x) = exp(-x)/(1 + exp(-x))^2   