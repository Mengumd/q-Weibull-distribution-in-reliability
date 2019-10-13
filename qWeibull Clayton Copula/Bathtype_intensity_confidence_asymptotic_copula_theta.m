global lifes;
global Cens;

[TBF,~,~] = xlsread('bathtype_intensity.xlsx','sheet7','A1:A44');
[Cens,~,~] = xlsread('bathtype_intensity.xlsx','sheet7','C1:C44');

lifes = zeros(1,length(TBF));
lifes(1) = TBF(1);
for i = 2:length(TBF)
    lifes(i) = lifes(i-1) + TBF(i) ;
end

[Obj, sol, output]=ABC_direct_throw_local_search_adaptive_copula_theta(1);


I = zeros(13,13);
ACI = zeros(1,26);

t = lifes;

syms theta loglik_vector Copula R
syms beta_1 beta_2 beta_3 beta_4 beta_5 beta_6
syms eta_1 eta_2 eta_3 eta_4 eta_5 eta_6


for ii = 1:length(t)
    
    
    R(1) = exp(-(t(ii)./eta_1).^beta_1) ;
    R(2) = exp(-(t(ii)./eta_2).^beta_2) ;
    R(3) = exp(-(t(ii)./eta_3).^beta_3) ;
    R(4) = exp(-(t(ii)./eta_4).^beta_4) ;
    R(5) = exp(-(t(ii)./eta_5).^beta_5) ;
    R(6) = exp(-(t(ii)./eta_6).^beta_6) ;
        
      
%     Copula(ii) = ( max(sum(R.^theta)-6+1, 0) )^(1/theta) ; 
    Copula(ii) = ( sum(R.^theta)-6+1 )^(1/theta) ; 
end


for ii = 1:length(t)       
    
    kk = Cens(ii);        
    
    if kk == 1
        loglik_vector(ii) = log(beta_1/(eta_1^beta_1)) + (beta_1-1)*log(t(ii)) - theta*(t(ii)/eta_1)^beta_1 - theta * log(Copula(ii));
    end
    
    if kk == 2
        loglik_vector(ii) = log(beta_2/(eta_2^beta_2)) + (beta_2-1)*log(t(ii)) - theta*(t(ii)/eta_2)^beta_2 - theta * log(Copula(ii));
    end
    
    if kk == 3
        loglik_vector(ii) = log(beta_3/(eta_3^beta_3)) + (beta_3-1)*log(t(ii)) - theta*(t(ii)/eta_3)^beta_3 - theta * log(Copula(ii));
    end
    
    if kk == 4
        loglik_vector(ii) = log(beta_4/(eta_4^beta_4)) + (beta_4-1)*log(t(ii)) - theta*(t(ii)/eta_4)^beta_4 - theta * log(Copula(ii));
    end
    
    if kk == 5
        loglik_vector(ii) = log(beta_5/(eta_5^beta_5)) + (beta_5-1)*log(t(ii)) - theta*(t(ii)/eta_5)^beta_5 - theta * log(Copula(ii));
    end
    
    if kk == 6
        loglik_vector(ii) = log(beta_6/(eta_6^beta_6)) + (beta_6-1)*log(t(ii)) - theta*(t(ii)/eta_6)^beta_6 - theta * log(Copula(ii));
    end
       
end

loglik = sum(loglik_vector) + log(Copula(end));


theta_loglik = diff(loglik, theta, 2);
beta_1_loglik = diff(loglik, beta_1, 2);
beta_2_loglik = diff(loglik, beta_2, 2);
beta_3_loglik = diff(loglik, beta_3, 2);
beta_4_loglik = diff(loglik, beta_4, 2);
beta_5_loglik = diff(loglik, beta_5, 2);
beta_6_loglik = diff(loglik, beta_6, 2);
eta_1_loglik = diff(loglik, eta_1, 2);
eta_2_loglik = diff(loglik, eta_2, 2);
eta_3_loglik = diff(loglik, eta_3, 2);
eta_4_loglik = diff(loglik, eta_4, 2);
eta_5_loglik = diff(loglik, eta_5, 2);
eta_6_loglik = diff(loglik, eta_6, 2);
theta_beta_1_loglik = diff(loglik, theta, beta_1);
theta_beta_2_loglik = diff(loglik, theta, beta_2);
theta_beta_3_loglik = diff(loglik, theta, beta_3);
theta_beta_4_loglik = diff(loglik, theta, beta_4);
theta_beta_5_loglik = diff(loglik, theta, beta_5);
theta_beta_6_loglik = diff(loglik, theta, beta_6);
theta_eta_1_loglik = diff(loglik, theta, eta_1);
theta_eta_2_loglik = diff(loglik, theta, eta_2);
theta_eta_3_loglik = diff(loglik, theta, eta_3);
theta_eta_4_loglik = diff(loglik, theta, eta_4);
theta_eta_5_loglik = diff(loglik, theta, eta_5);
theta_eta_6_loglik = diff(loglik, theta, eta_6);
beta_1_beta_2_loglik = diff(loglik, beta_1, beta_2);
beta_1_beta_3_loglik = diff(loglik, beta_1, beta_3);
beta_1_beta_4_loglik = diff(loglik, beta_1, beta_4);
beta_1_beta_5_loglik = diff(loglik, beta_1, beta_5);
beta_1_beta_6_loglik = diff(loglik, beta_1, beta_6);
beta_1_eta_1_loglik = diff(loglik, beta_1, eta_1);
beta_1_eta_2_loglik = diff(loglik, beta_1, eta_2);
beta_1_eta_3_loglik = diff(loglik, beta_1, eta_3);
beta_1_eta_4_loglik = diff(loglik, beta_1, eta_4);
beta_1_eta_5_loglik = diff(loglik, beta_1, eta_5);
beta_1_eta_6_loglik = diff(loglik, beta_1, eta_6);
beta_2_beta_3_loglik = diff(loglik, beta_2, beta_3);
beta_2_beta_4_loglik = diff(loglik, beta_2, beta_4);
beta_2_beta_5_loglik = diff(loglik, beta_2, beta_5);
beta_2_beta_6_loglik = diff(loglik, beta_2, beta_6);
beta_2_eta_1_loglik = diff(loglik, beta_2, eta_1);
beta_2_eta_2_loglik = diff(loglik, beta_2, eta_2);
beta_2_eta_3_loglik = diff(loglik, beta_2, eta_3);
beta_2_eta_4_loglik = diff(loglik, beta_2, eta_4);
beta_2_eta_5_loglik = diff(loglik, beta_2, eta_5);
beta_2_eta_6_loglik = diff(loglik, beta_2, eta_6);
beta_3_beta_4_loglik = diff(loglik, beta_3, beta_4);
beta_3_beta_5_loglik = diff(loglik, beta_3, beta_5);
beta_3_beta_6_loglik = diff(loglik, beta_3, beta_6);
beta_3_eta_1_loglik = diff(loglik, beta_3, eta_1);
beta_3_eta_2_loglik = diff(loglik, beta_3, eta_2);
beta_3_eta_3_loglik = diff(loglik, beta_3, eta_3);
beta_3_eta_4_loglik = diff(loglik, beta_3, eta_4);
beta_3_eta_5_loglik = diff(loglik, beta_3, eta_5);
beta_3_eta_6_loglik = diff(loglik, beta_3, eta_6);
beta_4_beta_5_loglik = diff(loglik, beta_4, beta_5);
beta_4_beta_6_loglik = diff(loglik, beta_4, beta_6);
beta_4_eta_1_loglik = diff(loglik, beta_4, eta_1);
beta_4_eta_2_loglik = diff(loglik, beta_4, eta_2);
beta_4_eta_3_loglik = diff(loglik, beta_4, eta_3);
beta_4_eta_4_loglik = diff(loglik, beta_4, eta_4);
beta_4_eta_5_loglik = diff(loglik, beta_4, eta_5);
beta_4_eta_6_loglik = diff(loglik, beta_4, eta_6);
beta_5_beta_6_loglik = diff(loglik, beta_5, beta_6);
beta_5_eta_1_loglik = diff(loglik, beta_5, eta_1);
beta_5_eta_2_loglik = diff(loglik, beta_5, eta_2);
beta_5_eta_3_loglik = diff(loglik, beta_5, eta_3);
beta_5_eta_4_loglik = diff(loglik, beta_5, eta_4);
beta_5_eta_5_loglik = diff(loglik, beta_5, eta_5);
beta_5_eta_6_loglik = diff(loglik, beta_5, eta_6);
beta_6_eta_1_loglik = diff(loglik, beta_6, eta_1);
beta_6_eta_2_loglik = diff(loglik, beta_6, eta_2);
beta_6_eta_3_loglik = diff(loglik, beta_6, eta_3);
beta_6_eta_4_loglik = diff(loglik, beta_6, eta_4);
beta_6_eta_5_loglik = diff(loglik, beta_6, eta_5);
beta_6_eta_6_loglik = diff(loglik, beta_6, eta_6);
eta_1_eta_2_loglik = diff(loglik, eta_1, eta_2);
eta_1_eta_3_loglik = diff(loglik, eta_1, eta_3);
eta_1_eta_4_loglik = diff(loglik, eta_1, eta_4);
eta_1_eta_5_loglik = diff(loglik, eta_1, eta_5);
eta_1_eta_6_loglik = diff(loglik, eta_1, eta_6);
eta_2_eta_3_loglik = diff(loglik, eta_2, eta_3);
eta_2_eta_4_loglik = diff(loglik, eta_2, eta_4);
eta_2_eta_5_loglik = diff(loglik, eta_2, eta_5);
eta_2_eta_6_loglik = diff(loglik, eta_2, eta_6);
eta_3_eta_4_loglik = diff(loglik, eta_3, eta_4);
eta_3_eta_5_loglik = diff(loglik, eta_3, eta_5);
eta_3_eta_6_loglik = diff(loglik, eta_3, eta_6);
eta_4_eta_5_loglik = diff(loglik, eta_4, eta_5);
eta_4_eta_6_loglik = diff(loglik, eta_4, eta_6);
eta_5_eta_6_loglik = diff(loglik, eta_5, eta_6);


theta = sol(1);
beta_1 = sol(1,2);
beta_2 = sol(1,3);
beta_3 = sol(1,4);
beta_4 = sol(1,5);
beta_5 = sol(1,6);
beta_6 = sol(1,7);
eta_1 = sol(1,8);
eta_2 = sol(1,9);
eta_3 = sol(1,10);
eta_4 = sol(1,11);
eta_5 = sol(1,12);
eta_6 = sol(1,13);


theta_loglik = subs(theta_loglik);
beta_1_loglik = subs(beta_1_loglik);
beta_2_loglik = subs(beta_2_loglik);
beta_3_loglik = subs(beta_3_loglik);
beta_4_loglik = subs(beta_4_loglik);
beta_5_loglik = subs(beta_5_loglik);
beta_6_loglik = subs(beta_6_loglik);
eta_1_loglik = subs(eta_1_loglik);
eta_2_loglik = subs(eta_2_loglik);
eta_3_loglik = subs(eta_3_loglik);
eta_4_loglik = subs(eta_4_loglik);
eta_5_loglik = subs(eta_5_loglik);
eta_6_loglik = subs(eta_6_loglik);
theta_beta_1_loglik = subs(theta_beta_1_loglik);
theta_beta_2_loglik = subs(theta_beta_2_loglik);
theta_beta_3_loglik = subs(theta_beta_3_loglik);
theta_beta_4_loglik = subs(theta_beta_4_loglik);
theta_beta_5_loglik = subs(theta_beta_5_loglik);
theta_beta_6_loglik = subs(theta_beta_6_loglik);
theta_eta_1_loglik = subs(theta_eta_1_loglik);
theta_eta_2_loglik = subs(theta_eta_2_loglik);
theta_eta_3_loglik = subs(theta_eta_3_loglik);
theta_eta_4_loglik = subs(theta_eta_4_loglik);
theta_eta_5_loglik = subs(theta_eta_5_loglik);
theta_eta_6_loglik = subs(theta_eta_6_loglik);
beta_1_beta_2_loglik = subs(beta_1_beta_2_loglik);
beta_1_beta_3_loglik = subs(beta_1_beta_3_loglik);
beta_1_beta_4_loglik = subs(beta_1_beta_4_loglik);
beta_1_beta_5_loglik = subs(beta_1_beta_5_loglik);
beta_1_beta_6_loglik = subs(beta_1_beta_6_loglik);
beta_1_eta_1_loglik = subs(beta_1_eta_1_loglik);
beta_1_eta_2_loglik = subs(beta_1_eta_2_loglik);
beta_1_eta_3_loglik = subs(beta_1_eta_3_loglik);
beta_1_eta_4_loglik = subs(beta_1_eta_4_loglik);
beta_1_eta_5_loglik = subs(beta_1_eta_5_loglik);
beta_1_eta_6_loglik = subs(beta_1_eta_6_loglik);
beta_2_beta_3_loglik = subs(beta_2_beta_3_loglik);
beta_2_beta_4_loglik = subs(beta_2_beta_4_loglik);
beta_2_beta_5_loglik = subs(beta_2_beta_5_loglik);
beta_2_beta_6_loglik = subs(beta_2_beta_6_loglik);
beta_2_eta_1_loglik = subs(beta_2_eta_1_loglik);
beta_2_eta_2_loglik = subs(beta_2_eta_2_loglik);
beta_2_eta_3_loglik = subs(beta_2_eta_3_loglik);
beta_2_eta_4_loglik = subs(beta_2_eta_4_loglik);
beta_2_eta_5_loglik = subs(beta_2_eta_5_loglik);
beta_2_eta_6_loglik = subs(beta_2_eta_6_loglik);
beta_3_beta_4_loglik = subs(beta_3_beta_4_loglik);
beta_3_beta_5_loglik = subs(beta_3_beta_5_loglik);
beta_3_beta_6_loglik = subs(beta_3_beta_6_loglik);
beta_3_eta_1_loglik = subs(beta_3_eta_1_loglik);
beta_3_eta_2_loglik = subs(beta_3_eta_2_loglik);
beta_3_eta_3_loglik = subs(beta_3_eta_3_loglik);
beta_3_eta_4_loglik = subs(beta_3_eta_4_loglik);
beta_3_eta_5_loglik = subs(beta_3_eta_5_loglik);
beta_3_eta_6_loglik = subs(beta_3_eta_6_loglik);
beta_4_beta_5_loglik = subs(beta_4_beta_5_loglik);
beta_4_beta_6_loglik = subs(beta_4_beta_6_loglik);
beta_4_eta_1_loglik = subs(beta_4_eta_1_loglik);
beta_4_eta_2_loglik = subs(beta_4_eta_2_loglik);
beta_4_eta_3_loglik = subs(beta_4_eta_3_loglik);
beta_4_eta_4_loglik = subs(beta_4_eta_4_loglik);
beta_4_eta_5_loglik = subs(beta_4_eta_5_loglik);
beta_4_eta_6_loglik = subs(beta_4_eta_6_loglik);
beta_5_beta_6_loglik = subs(beta_5_beta_6_loglik);
beta_5_eta_1_loglik = subs(beta_5_eta_1_loglik);
beta_5_eta_2_loglik = subs(beta_5_eta_2_loglik);
beta_5_eta_3_loglik = subs(beta_5_eta_3_loglik);
beta_5_eta_4_loglik = subs(beta_5_eta_4_loglik);
beta_5_eta_5_loglik = subs(beta_5_eta_5_loglik);
beta_5_eta_6_loglik = subs(beta_5_eta_6_loglik);
beta_6_eta_1_loglik = subs(beta_6_eta_1_loglik);
beta_6_eta_2_loglik = subs(beta_6_eta_2_loglik);
beta_6_eta_3_loglik = subs(beta_6_eta_3_loglik);
beta_6_eta_4_loglik = subs(beta_6_eta_4_loglik);
beta_6_eta_5_loglik = subs(beta_6_eta_5_loglik);
beta_6_eta_6_loglik = subs(beta_6_eta_6_loglik);
eta_1_eta_2_loglik = subs(eta_1_eta_2_loglik);
eta_1_eta_3_loglik = subs(eta_1_eta_3_loglik);
eta_1_eta_4_loglik = subs(eta_1_eta_4_loglik);
eta_1_eta_5_loglik = subs(eta_1_eta_5_loglik);
eta_1_eta_6_loglik = subs(eta_1_eta_6_loglik);
eta_2_eta_3_loglik = subs(eta_2_eta_3_loglik);
eta_2_eta_4_loglik = subs(eta_2_eta_4_loglik);
eta_2_eta_5_loglik = subs(eta_2_eta_5_loglik);
eta_2_eta_6_loglik = subs(eta_2_eta_6_loglik);
eta_3_eta_4_loglik = subs(eta_3_eta_4_loglik);
eta_3_eta_5_loglik = subs(eta_3_eta_5_loglik);
eta_3_eta_6_loglik = subs(eta_3_eta_6_loglik);
eta_4_eta_5_loglik = subs(eta_4_eta_5_loglik);
eta_4_eta_6_loglik = subs(eta_4_eta_6_loglik);
eta_5_eta_6_loglik = subs(eta_5_eta_6_loglik);



theta_loglik = double(theta_loglik);
beta_1_loglik = double(beta_1_loglik);
beta_2_loglik = double(beta_2_loglik);
beta_3_loglik = double(beta_3_loglik);
beta_4_loglik = double(beta_4_loglik);
beta_5_loglik = double(beta_5_loglik);
beta_6_loglik = double(beta_6_loglik);
eta_1_loglik = double(eta_1_loglik);
eta_2_loglik = double(eta_2_loglik);
eta_3_loglik = double(eta_3_loglik);
eta_4_loglik = double(eta_4_loglik);
eta_5_loglik = double(eta_5_loglik);
eta_6_loglik = double(eta_6_loglik);
theta_beta_1_loglik = double(theta_beta_1_loglik);
theta_beta_2_loglik = double(theta_beta_2_loglik);
theta_beta_3_loglik = double(theta_beta_3_loglik);
theta_beta_4_loglik = double(theta_beta_4_loglik);
theta_beta_5_loglik = double(theta_beta_5_loglik);
theta_beta_6_loglik = double(theta_beta_6_loglik);
theta_eta_1_loglik = double(theta_eta_1_loglik);
theta_eta_2_loglik = double(theta_eta_2_loglik);
theta_eta_3_loglik = double(theta_eta_3_loglik);
theta_eta_4_loglik = double(theta_eta_4_loglik);
theta_eta_5_loglik = double(theta_eta_5_loglik);
theta_eta_6_loglik = double(theta_eta_6_loglik);
beta_1_beta_2_loglik = double(beta_1_beta_2_loglik);
beta_1_beta_3_loglik = double(beta_1_beta_3_loglik);
beta_1_beta_4_loglik = double(beta_1_beta_4_loglik);
beta_1_beta_5_loglik = double(beta_1_beta_5_loglik);
beta_1_beta_6_loglik = double(beta_1_beta_6_loglik);
beta_1_eta_1_loglik = double(beta_1_eta_1_loglik);
beta_1_eta_2_loglik = double(beta_1_eta_2_loglik);
beta_1_eta_3_loglik = double(beta_1_eta_3_loglik);
beta_1_eta_4_loglik = double(beta_1_eta_4_loglik);
beta_1_eta_5_loglik = double(beta_1_eta_5_loglik);
beta_1_eta_6_loglik = double(beta_1_eta_6_loglik);
beta_2_beta_3_loglik = double(beta_2_beta_3_loglik);
beta_2_beta_4_loglik = double(beta_2_beta_4_loglik);
beta_2_beta_5_loglik = double(beta_2_beta_5_loglik);
beta_2_beta_6_loglik = double(beta_2_beta_6_loglik);
beta_2_eta_1_loglik = double(beta_2_eta_1_loglik);
beta_2_eta_2_loglik = double(beta_2_eta_2_loglik);
beta_2_eta_3_loglik = double(beta_2_eta_3_loglik);
beta_2_eta_4_loglik = double(beta_2_eta_4_loglik);
beta_2_eta_5_loglik = double(beta_2_eta_5_loglik);
beta_2_eta_6_loglik = double(beta_2_eta_6_loglik);
beta_3_beta_4_loglik = double(beta_3_beta_4_loglik);
beta_3_beta_5_loglik = double(beta_3_beta_5_loglik);
beta_3_beta_6_loglik = double(beta_3_beta_6_loglik);
beta_3_eta_1_loglik = double(beta_3_eta_1_loglik);
beta_3_eta_2_loglik = double(beta_3_eta_2_loglik);
beta_3_eta_3_loglik = double(beta_3_eta_3_loglik);
beta_3_eta_4_loglik = double(beta_3_eta_4_loglik);
beta_3_eta_5_loglik = double(beta_3_eta_5_loglik);
beta_3_eta_6_loglik = double(beta_3_eta_6_loglik);
beta_4_beta_5_loglik = double(beta_4_beta_5_loglik);
beta_4_beta_6_loglik = double(beta_4_beta_6_loglik);
beta_4_eta_1_loglik = double(beta_4_eta_1_loglik);
beta_4_eta_2_loglik = double(beta_4_eta_2_loglik);
beta_4_eta_3_loglik = double(beta_4_eta_3_loglik);
beta_4_eta_4_loglik = double(beta_4_eta_4_loglik);
beta_4_eta_5_loglik = double(beta_4_eta_5_loglik);
beta_4_eta_6_loglik = double(beta_4_eta_6_loglik);
beta_5_beta_6_loglik = double(beta_5_beta_6_loglik);
beta_5_eta_1_loglik = double(beta_5_eta_1_loglik);
beta_5_eta_2_loglik = double(beta_5_eta_2_loglik);
beta_5_eta_3_loglik = double(beta_5_eta_3_loglik);
beta_5_eta_4_loglik = double(beta_5_eta_4_loglik);
beta_5_eta_5_loglik = double(beta_5_eta_5_loglik);
beta_5_eta_6_loglik = double(beta_5_eta_6_loglik);
beta_6_eta_1_loglik = double(beta_6_eta_1_loglik);
beta_6_eta_2_loglik = double(beta_6_eta_2_loglik);
beta_6_eta_3_loglik = double(beta_6_eta_3_loglik);
beta_6_eta_4_loglik = double(beta_6_eta_4_loglik);
beta_6_eta_5_loglik = double(beta_6_eta_5_loglik);
beta_6_eta_6_loglik = double(beta_6_eta_6_loglik);
eta_1_eta_2_loglik = double(eta_1_eta_2_loglik);
eta_1_eta_3_loglik = double(eta_1_eta_3_loglik);
eta_1_eta_4_loglik = double(eta_1_eta_4_loglik);
eta_1_eta_5_loglik = double(eta_1_eta_5_loglik);
eta_1_eta_6_loglik = double(eta_1_eta_6_loglik);
eta_2_eta_3_loglik = double(eta_2_eta_3_loglik);
eta_2_eta_4_loglik = double(eta_2_eta_4_loglik);
eta_2_eta_5_loglik = double(eta_2_eta_5_loglik);
eta_2_eta_6_loglik = double(eta_2_eta_6_loglik);
eta_3_eta_4_loglik = double(eta_3_eta_4_loglik);
eta_3_eta_5_loglik = double(eta_3_eta_5_loglik);
eta_3_eta_6_loglik = double(eta_3_eta_6_loglik);
eta_4_eta_5_loglik = double(eta_4_eta_5_loglik);
eta_4_eta_6_loglik = double(eta_4_eta_6_loglik);
eta_5_eta_6_loglik = double(eta_5_eta_6_loglik);




I(1,1) = theta_loglik;
I(2,2) = beta_1_loglik;
I(3,3) = beta_2_loglik;
I(4,4) = beta_3_loglik;
I(5,5) = beta_4_loglik;
I(6,6) = beta_5_loglik;
I(7,7) = beta_6_loglik;
I(8,8) = eta_1_loglik;
I(9,9) = eta_2_loglik;
I(10,10) = eta_3_loglik;
I(11,11) = eta_4_loglik;
I(12,12) = eta_5_loglik;
I(13,13) = eta_6_loglik;
I(1,2) = theta_beta_1_loglik;
I(2,1) = theta_beta_1_loglik;
I(1,3) = theta_beta_2_loglik;
I(3,1) = theta_beta_2_loglik;
I(1,4) = theta_beta_3_loglik;
I(4,1) = theta_beta_3_loglik;
I(1,5) = theta_beta_4_loglik;
I(5,1) = theta_beta_4_loglik;
I(1,6) = theta_beta_5_loglik;
I(6,1) = theta_beta_5_loglik;
I(1,7) = theta_beta_6_loglik;
I(7,1) = theta_beta_6_loglik;
I(1,8) = theta_eta_1_loglik;
I(8,1) = theta_eta_1_loglik;
I(1,9) = theta_eta_2_loglik;
I(9,1) = theta_eta_2_loglik;
I(1,10) = theta_eta_3_loglik;
I(10,1) = theta_eta_3_loglik;
I(1,11) = theta_eta_4_loglik;
I(11,1) = theta_eta_4_loglik;
I(1,12) = theta_eta_5_loglik;
I(12,1) = theta_eta_5_loglik;
I(1,13) = theta_eta_6_loglik;
I(13,1) = theta_eta_6_loglik;
I(2,3) = beta_1_beta_2_loglik;
I(3,2) = beta_1_beta_2_loglik;
I(2,4) = beta_1_beta_3_loglik;
I(4,2) = beta_1_beta_3_loglik;
I(2,5) = beta_1_beta_4_loglik;
I(5,2) = beta_1_beta_4_loglik;
I(2,6) = beta_1_beta_5_loglik;
I(6,2) = beta_1_beta_5_loglik;
I(2,7) = beta_1_beta_6_loglik;
I(7,2) = beta_1_beta_6_loglik;
I(2,8) = beta_1_eta_1_loglik;
I(8,2) = beta_1_eta_1_loglik;
I(2,9) = beta_1_eta_2_loglik;
I(9,2) = beta_1_eta_2_loglik;
I(2,10) = beta_1_eta_3_loglik;
I(10,2) = beta_1_eta_3_loglik;
I(2,11) = beta_1_eta_4_loglik;
I(11,2) = beta_1_eta_4_loglik;
I(2,12) = beta_1_eta_5_loglik;
I(12,2) = beta_1_eta_5_loglik;
I(2,13) = beta_1_eta_6_loglik;
I(13,2) = beta_1_eta_6_loglik;
I(3,4) = beta_2_beta_3_loglik;
I(4,3) = beta_2_beta_3_loglik;
I(3,5) = beta_2_beta_4_loglik;
I(5,3) = beta_2_beta_4_loglik;
I(3,6) = beta_2_beta_5_loglik;
I(6,3) = beta_2_beta_5_loglik;
I(3,7) = beta_2_beta_6_loglik;
I(7,3) = beta_2_beta_6_loglik;
I(3,8) = beta_2_eta_1_loglik;
I(8,3) = beta_2_eta_1_loglik;
I(3,9) = beta_2_eta_2_loglik;
I(9,3) = beta_2_eta_2_loglik;
I(3,10) = beta_2_eta_3_loglik;
I(10,3) = beta_2_eta_3_loglik;
I(3,11) = beta_2_eta_4_loglik;
I(11,3) = beta_2_eta_4_loglik;
I(3,12) = beta_2_eta_5_loglik;
I(12,3) = beta_2_eta_5_loglik;
I(3,13) = beta_2_eta_6_loglik;
I(13,3) = beta_2_eta_6_loglik;
I(4,5) = beta_3_beta_4_loglik;
I(5,4) = beta_3_beta_4_loglik;
I(4,6) = beta_3_beta_5_loglik;
I(6,4) = beta_3_beta_5_loglik;
I(4,7) = beta_3_beta_6_loglik;
I(7,4) = beta_3_beta_6_loglik;
I(4,8) = beta_3_eta_1_loglik;
I(8,4) = beta_3_eta_1_loglik;
I(4,9) = beta_3_eta_2_loglik;
I(9,4) = beta_3_eta_2_loglik;
I(4,10) = beta_3_eta_3_loglik;
I(10,4) = beta_3_eta_3_loglik;
I(4,11) = beta_3_eta_4_loglik;
I(11,4) = beta_3_eta_4_loglik;
I(4,12) = beta_3_eta_5_loglik;
I(12,4) = beta_3_eta_5_loglik;
I(4,13) = beta_3_eta_6_loglik;
I(13,4) = beta_3_eta_6_loglik;
I(5,6) = beta_4_beta_5_loglik;
I(6,5) = beta_4_beta_5_loglik;
I(5,7) = beta_4_beta_6_loglik;
I(7,5) = beta_4_beta_6_loglik;
I(5,8) = beta_4_eta_1_loglik;
I(8,5) = beta_4_eta_1_loglik;
I(5,9) = beta_4_eta_2_loglik;
I(9,5) = beta_4_eta_2_loglik;
I(5,10) = beta_4_eta_3_loglik;
I(10,5) = beta_4_eta_3_loglik;
I(5,11) = beta_4_eta_4_loglik;
I(11,5) = beta_4_eta_4_loglik;
I(5,12) = beta_4_eta_5_loglik;
I(12,5) = beta_4_eta_5_loglik;
I(5,13) = beta_4_eta_6_loglik;
I(13,5) = beta_4_eta_6_loglik;
I(6,7) = beta_5_beta_6_loglik;
I(7,6) = beta_5_beta_6_loglik;
I(6,8) = beta_5_eta_1_loglik;
I(8,6) = beta_5_eta_1_loglik;
I(6,9) = beta_5_eta_2_loglik;
I(9,6) = beta_5_eta_2_loglik;
I(6,10) = beta_5_eta_3_loglik;
I(10,6) = beta_5_eta_3_loglik;
I(6,11) = beta_5_eta_4_loglik;
I(11,6) = beta_5_eta_4_loglik;
I(6,12) = beta_5_eta_5_loglik;
I(12,6) = beta_5_eta_5_loglik;
I(6,13) = beta_5_eta_6_loglik;
I(13,6) = beta_5_eta_6_loglik;
I(7,8) = beta_6_eta_1_loglik;
I(8,7) = beta_6_eta_1_loglik;
I(7,9) = beta_6_eta_2_loglik;
I(9,7) = beta_6_eta_2_loglik;
I(7,10) = beta_6_eta_3_loglik;
I(10,7) = beta_6_eta_3_loglik;
I(7,11) = beta_6_eta_4_loglik;
I(11,7) = beta_6_eta_4_loglik;
I(7,12) = beta_6_eta_5_loglik;
I(12,7) = beta_6_eta_5_loglik;
I(7,13) = beta_6_eta_6_loglik;
I(13,7) = beta_6_eta_6_loglik;
I(8,9) = eta_1_eta_2_loglik;
I(9,8) = eta_1_eta_2_loglik;
I(8,10) = eta_1_eta_3_loglik;
I(10,8) = eta_1_eta_3_loglik;
I(8,11) = eta_1_eta_4_loglik;
I(11,8) = eta_1_eta_4_loglik;
I(8,12) = eta_1_eta_5_loglik;
I(12,8) = eta_1_eta_5_loglik;
I(8,13) = eta_1_eta_6_loglik;
I(13,8) = eta_1_eta_6_loglik;
I(9,10) = eta_2_eta_3_loglik;
I(10,9) = eta_2_eta_3_loglik;
I(9,11) = eta_2_eta_4_loglik;
I(11,9) = eta_2_eta_4_loglik;
I(9,12) = eta_2_eta_5_loglik;
I(12,9) = eta_2_eta_5_loglik;
I(9,13) = eta_2_eta_6_loglik;
I(13,9) = eta_2_eta_6_loglik;
I(10,11) = eta_3_eta_4_loglik;
I(11,10) = eta_3_eta_4_loglik;
I(10,12) = eta_3_eta_5_loglik;
I(12,10) = eta_3_eta_5_loglik;
I(10,13) = eta_3_eta_6_loglik;
I(13,10) = eta_3_eta_6_loglik;
I(11,12) = eta_4_eta_5_loglik;
I(12,11) = eta_4_eta_5_loglik;
I(11,13) = eta_4_eta_6_loglik;
I(13,11) = eta_4_eta_6_loglik;
I(12,13) = eta_5_eta_6_loglik;
I(13,12) = eta_5_eta_6_loglik;



var = -I^(-1);
var_theta = var(1,1);
var_beta_1 = var(2,2);
var_beta_2 = var(3,3);
var_beta_3 = var(4,4);
var_beta_4 = var(5,5);
var_beta_5 = var(6,6);
var_beta_6 = var(7,7);
var_eta_1 = var(8,8);
var_eta_2 = var(9,9);
var_eta_3 = var(10,10);
var_eta_4 = var(11,11);
var_eta_5 = var(12,12);
var_eta_6 = var(13,13);

Z1 = -1.64;
Z2 = 1.64;

theta_ub = sol(1) + Z2*(var_theta)^(0.5);
theta_lb = sol(1) + Z1*(var_theta)^(0.5);
beta_1_ub = sol(1,2) + Z2*(var_beta_1)^(0.5);
beta_1_lb = sol(1,2) + Z1*(var_beta_1)^(0.5);
beta_2_ub = sol(1,3) + Z2*(var_beta_2)^(0.5);
beta_2_lb = sol(1,3) + Z1*(var_beta_2)^(0.5);
beta_3_ub = sol(1,4) + Z2*(var_beta_3)^(0.5);
beta_3_lb = sol(1,4) + Z1*(var_beta_3)^(0.5);
beta_4_ub = sol(1,5) + Z2*(var_beta_4)^(0.5);
beta_4_lb = sol(1,5) + Z1*(var_beta_4)^(0.5);
beta_5_ub = sol(1,6) + Z2*(var_beta_5)^(0.5);
beta_5_lb = sol(1,6) + Z1*(var_beta_5)^(0.5);
beta_6_ub = sol(1,7) + Z2*(var_beta_6)^(0.5);
beta_6_lb = sol(1,7) + Z1*(var_beta_6)^(0.5);
eta_1_ub = sol(1,8) + Z2*(var_eta_1)^(0.5);
eta_1_lb = sol(1,8) + Z1*(var_eta_1)^(0.5);
eta_2_ub = sol(1,9) + Z2*(var_eta_2)^(0.5);
eta_2_lb = sol(1,9) + Z1*(var_eta_2)^(0.5);
eta_3_ub = sol(1,10) + Z2*(var_eta_3)^(0.5);
eta_3_lb = sol(1,10) + Z1*(var_eta_3)^(0.5);
eta_4_ub = sol(1,11) + Z2*(var_eta_4)^(0.5);
eta_4_lb = sol(1,11) + Z1*(var_eta_4)^(0.5);
eta_5_ub = sol(1,12) + Z2*(var_eta_5)^(0.5);
eta_5_lb = sol(1,12) + Z1*(var_eta_5)^(0.5);
eta_6_ub = sol(1,13) + Z2*(var_eta_6)^(0.5);
eta_6_lb = sol(1,13) + Z1*(var_eta_6)^(0.5);


theta_ub = sol(1) + Z2*(var_theta)^(0.5);
theta_lb = sol(1) + Z1*(var_theta)^(0.5);
ACI(1,1) = beta_1_ub;
ACI(2,1) = beta_1_lb;
ACI(1,2) = beta_2_ub;
ACI(2,2) = beta_2_lb;
ACI(1,3) = beta_3_ub;
ACI(2,3) = beta_3_lb;
ACI(1,4) = beta_4_ub;
ACI(2,4) = beta_4_lb;
ACI(1,5) = beta_5_ub;
ACI(2,5) = beta_5_lb;
ACI(1,6) = beta_6_ub;
ACI(2,6) = beta_6_lb;
ACI(1,7) = eta_1_ub;
ACI(2,7) = eta_1_lb;
ACI(1,8) = eta_2_ub;
ACI(2,8) = eta_2_lb;
ACI(1,9) = eta_3_ub;
ACI(2,9) = eta_3_lb;
ACI(1,10) = eta_4_ub;
ACI(2,10) = eta_4_lb;
ACI(1,11) = eta_5_ub;
ACI(2,11) = eta_5_lb;
ACI(1,12) = eta_6_ub;
ACI(2,12) = eta_6_lb;




   
    
if Copula(end) == 0
    loglik = -inf;    
    return
end

