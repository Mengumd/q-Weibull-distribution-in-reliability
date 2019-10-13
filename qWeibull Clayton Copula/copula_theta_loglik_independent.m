function loglik = copula_theta_loglik_independent(paras)
%% this function is to calculate the log of likelyhood, with a given paras
% paras: the parameters,i.e. [q, beta, yita]
% lifes: experimental test life data
% global Para_beta;
% global Para_eta;

Para_beta = paras(1:6);
Para_eta = paras(7:12);

global lifes;
t = lifes;

global Cens;

loglik_vector = zeros(1,length(lifes));
      
R = exp(-(t(end)./Para_eta).^Para_beta) ;    
Copula = prod(R) ;

if Copula == 0
    loglik = -inf;    
    return
end


for ii = 1:length(t)       
    if Cens(ii) ~= 0
        kk = Cens(ii);
        loglik_vector(ii) = log(Para_beta(kk)/(Para_eta(kk)^Para_beta(kk))) + (Para_beta(kk)-1)*log(t(ii));       
    end        
end

loglik = sum(loglik_vector) + log(Copula);


end
