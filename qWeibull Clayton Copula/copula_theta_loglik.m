function loglik = copula_theta_loglik(paras)
%% this function is to calculate the log of likelyhood, with a given paras
% paras: the parameters,i.e. [q, beta, yita]
% lifes: experimental test life data
theta = paras(1);
% global Para_beta;
% global Para_eta;

Para_beta = paras(2:7);
Para_eta = paras(8:13);
d = length(Para_beta);

global lifes;
t = lifes;

global Cens;

Copula = zeros(1,length(lifes));
loglik_vector = zeros(1,length(lifes));

for ii = 1:length(t)       
    R = exp(-(t(ii)./Para_eta).^Para_beta) ;    
    Copula(ii) = ( max(sum(R.^theta)-d+1, 0) )^(1/theta) ;    
end


if Copula(end) == 0
    loglik = -inf;    
    return
end

for ii = 1:length(t)       
    if Cens(ii) ~= 0
        kk = Cens(ii);
        loglik_vector(ii) = log(Para_beta(kk)/(Para_eta(kk)^Para_beta(kk))) + (Para_beta(kk)-1)*log(t(ii)) - theta*(t(ii)/Para_eta(kk))^Para_beta(kk) - theta * log(Copula(ii));       
    end        
end

loglik = sum(loglik_vector) + log(Copula(end));


end
