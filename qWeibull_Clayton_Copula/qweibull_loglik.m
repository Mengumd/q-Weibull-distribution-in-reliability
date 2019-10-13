function loglik = qweibull_loglik(paras)
%% this function is to calculate the log of likelyhood, with a given paras
% paras: the parameters,i.e. [q, beta, yita]
% lifes: experimental test life data

q = paras(1);
beta = paras(2);
eta = paras(3);


global lifes;
t = lifes;

global Cens;

for ii = 1:length(t)
    if Cens(ii) == 0
        loglik_vector(ii) = log((2-q)*beta/(eta^beta)) + (beta-1)*log(t(ii)) - log(1-(1-q)*(t(ii)/eta).^beta) ;
    end
end

loglik = sum(loglik_vector) + ((2-q)/(1-q))*log(1-(1-q)*(t(end)/eta)^beta);
end


