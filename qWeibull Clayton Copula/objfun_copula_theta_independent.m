function y = objfun_copula_theta_independent(x)
global nFcn;
nFcn = nFcn + 1;
y = -copula_theta_loglik_independent(x);