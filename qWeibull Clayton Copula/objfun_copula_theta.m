function y = objfun_copula_theta(x)
global nFcn;
nFcn = nFcn + 1;
y = -copula_theta_loglik(x);