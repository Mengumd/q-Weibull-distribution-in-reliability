function loglik = copula_theta_loglik_hr(paras)
%% this function is to calculate the log of likelyhood, with a given paras
% paras: the parameters,i.e. [q, beta, yita]
% lifes: experimental test life data
theta = paras(1);
% global Para_beta;
% global Para_eta;

beta = paras(2:7);
eta = paras(8:13);
% beta = [2 2 0.5 0.5];
% eta = [5 5 5 5];
% beta = paras(2:7);
% eta = paras(8:13);
% beta = [0.77 0.44 3.66 0.86 1.19 2.32];
% eta = [114.7 20.3 1937 770 688 1914];
% beta = [2 2 2 2];
% eta = [5 5 5 5];
d = length(beta);


global lifes;
t = lifes;
global Cens;

inslog = zeros(1,length(lifes));
Copula = zeros(1,length(lifes));
loglik_vector = zeros(1,length(lifes));

for ii = 1:length(t)
    R = exp(-(t(ii)./eta).^beta) ;
    h = (beta./eta).*(t(ii)./eta).^(beta-1);
    inslog(ii) = sum((R.^theta).*h);
    Copula(ii) = ( max(sum(R.^theta)-d+1, 0) )^(1/theta) ;
end

if Copula(end) == 0
    loglik = -inf;  
    return
end


for ii = 1:length(t)
    if Cens(ii) == 0
        loglik_vector(ii) = log(inslog(ii))  - theta * log(Copula(ii));  
    end
end

loglik = sum(loglik_vector) + log(Copula(end));

if loglik == inf
    loglik = -inf;  
    return
end


end
