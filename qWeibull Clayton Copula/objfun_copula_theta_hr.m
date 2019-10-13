function y = objfun_copula_theta_hr(x)
global nFcn;
nFcn = nFcn + 1;
y = -copula_theta_loglik_hr(x);