function loglik = SPLP_loglik(paras)
%% this function is to calculate the log of likelyhood, with a given paras
% paras: the parameters,i.e. [q, beta, yita]
% lifes: experimental test life data
alpha_1 = paras(1);
beta_1 = paras(2);
alpha_2 = paras(3);
beta_2 = paras(4);

global lifes;
global Cens;

t = lifes;
for ii = 1:length(t)
    if Cens(ii) == 0
        loglik_vector(ii) = log(beta_1*t(ii).^(beta_1-1)/alpha_1^beta_1 + beta_2*t(ii).^(beta_2-1)/alpha_2^beta_2);
    end
end

loglik = sum(loglik_vector) - ((t(end)/alpha_1).^beta_1 + (t(end)/alpha_2).^beta_2);












end
