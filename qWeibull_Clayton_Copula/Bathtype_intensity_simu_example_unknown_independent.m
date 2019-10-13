
global lifes;
global Cens;

% [lifes,~,~] = xlsread('bathtype_intensity.xlsx','sheet2');
% 
% lifes = reshape(lifes',1,[]);
% lifes = lifes(~isnan(lifes));

[TBF,~,~] = xlsread('bathtype_intensity.xlsx','Sheet7','A1:A44');
% [Cens,~,~] = xlsread('bathtype_intensity.xlsx','Sheet7','B1:B44');

lifes = zeros(1,length(TBF));
lifes(1) = TBF(1);
for i = 2:length(TBF)
    lifes(i) = lifes(i-1) + TBF(i) ;
end

Cens = zeros(1,length(lifes));
% [Obj, sol, output]=ABC_direct_throw_local_search_adaptive_weibull(1);
[Obj, sol, output]=copula_theta_unknown_independent(1);
objMean = mean(Obj);
objStd = std(Obj);
solMean = mean(sol);
solStd = std(sol);


save Weibull_NHPP_sub_2;

q = solMean(1);
beta = solMean(2);
eta = solMean(3);
t = 0:1:2317;
% y = (2-q)*(beta/eta)*(t./eta).^(beta-1).*exp_q(-(t./eta).^beta,q);
% R = (1-(1-q)*(t./eta).^beta).^((2-q)/(1-q));
% h = y./R;

M_splp = ((t/alpha_1).^beta_1) + (t/alpha_2).^beta_2 ;

M_qWeibull = -((2-q)*log(1-(1-q)*(t/eta).^beta))/(1-q) ;

n = length(lifes);
Fn = zeros(1,n+1);
for j = 0:n
    Fn(j+1) = j  ;
end
M_data = [0 lifes];

figure;
plot(t,M_splp,'r--',t,M_qWeibull,'k-',M_data,Fn,'o');
legend('S-PLP','q-Weibull','Observed data');
title('Cumulative number of failures');
xlabel('Operating time t [in hours]');
ylabel('Cumulative number of failures');


axis([0 2200 0 45]);