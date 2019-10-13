function loglik = weibull_loglik(paras)
%% this function is to calculate the log of likelyhood, with a given paras
% paras: the parameters,i.e. [q, beta, yita]
% lifes: experimental test life data
beta = paras(1);
eta = paras(2);

global lifes;
t = lifes;

global Cens;

for ii = 1:length(t)
    if Cens(ii) == 0
        loglik_vector(ii) = log(beta/(eta^beta)) + (beta-1)*log(t(ii)) ;
    end
end

loglik = sum(loglik_vector) -(t(end)/eta)^beta;
end
