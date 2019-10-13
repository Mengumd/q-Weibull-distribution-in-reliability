
global lifes;

% [lifes,~,~] = xlsread('bathtype_intensity.xlsx','sheet2');
% 
% lifes = reshape(lifes',1,[]);
% lifes = lifes(~isnan(lifes));

[lifes,~,~] = xlsread('bathtype_intensity.xlsx','A2:A18');


[Obj, sol, output]=ABC_direct_throw_local_search_adaptive(1);
objMean = mean(Obj);
objStd = std(Obj);
solMean = mean(sol);
solStd = std(sol);


save qWeibull_NHPP_sub_1;

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