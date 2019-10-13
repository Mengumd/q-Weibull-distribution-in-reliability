function y = objfun_SPLP(x)
global nFcn;
nFcn = nFcn + 1;
y = -SPLP_loglik(x);