function y = objfun_weibull(x)
global nFcn;
nFcn = nFcn + 1;
y = -weibull_loglik(x);