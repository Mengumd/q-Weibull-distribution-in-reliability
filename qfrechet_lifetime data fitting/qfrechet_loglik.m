function loglik = qfrechet_loglik(paras)
%% this function is to calculate the log of likelyhood, with a given paras
% paras: the parameters,i.e. [q, beta, yita]
% lifes: experimental test life data
q = paras(1);  % q from pdf
kse = paras(2);
sigma = paras(3);

global lifes;
t = lifes;
loglik_vector = log(kse/sigma) + (-kse-1)*log(t/sigma) + log(exp_q((-1/(2-q))*(t/sigma).^(-kse),q));
loglik = sum(loglik_vector);
end
