function y = objfun_copula_theta_unknown_independent(x)
global nFcn;
nFcn = nFcn + 1;
y = -copula_theta_loglik_unknown_independent(x);