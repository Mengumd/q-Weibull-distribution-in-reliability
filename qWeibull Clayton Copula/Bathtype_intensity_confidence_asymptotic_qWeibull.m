global lifes;
global Cens;

% [lifes,~,~] = xlsread('bathtype_intensity.xlsx','sheet2');
% 
% lifes = reshape(lifes',1,[]);
% lifes = lifes(~isnan(lifes));

[TBF,~,~] = xlsread('bathtype_intensity.xlsx','sheet7','A1:A44');
[Cens,~,~] = xlsread('bathtype_intensity.xlsx','sheet7','B1:B44');

lifes = zeros(1,length(TBF));
lifes(1) = TBF(1);
for i = 2:length(TBF)
    lifes(i) = lifes(i-1) + TBF(i) ;
end

[Obj, sol, output]=ABC_direct_throw_local_search_adaptive(1);


I = zeros(3,3);
ACI = zeros(1,9);

t = lifes;

syms q beta eta loglik_vector

for ii = 1:length(t)
    if Cens(ii) == 0
        loglik_vector(ii) = log((2-q)*beta/(eta^beta)) + (beta-1)*log(t(ii)) - log(1-(1-q)*(t(ii)/eta).^beta) ;
    end
end

loglik = sum(loglik_vector) + ((2-q)/(1-q))*log(1-(1-q)*(t(end)/eta)^beta);

q_loglik = diff(loglik, q, 2);

beta_loglik = diff(loglik, beta, 2);

eta_loglik = diff(loglik, eta, 2);

q_beta_loglik = diff(loglik, q, beta);

beta_eta_loglik = diff(loglik, beta, eta);

eta_q_loglik = diff(loglik, q, eta);

q = sol(1);
beta = sol(1,2);
eta = sol(1,3);

q_loglik = subs(q_loglik);

beta_loglik = subs(beta_loglik);

eta_loglik = subs(eta_loglik);

q_beta_loglik = subs(q_beta_loglik);

beta_eta_loglik = subs(beta_eta_loglik);

eta_q_loglik = subs(eta_q_loglik);

I(1,1) = q_loglik;
I(2,2) = beta_loglik;
I(3,3) = eta_loglik;
I(1,2) = q_beta_loglik;
I(2,1) = q_beta_loglik;
I(1,3) = eta_q_loglik;
I(3,1) = eta_q_loglik;
I(2,3) = beta_eta_loglik;
I(3,2) = beta_eta_loglik;

var = -I^(-1);
var_q = var(1,1);
var_beta = var(2,2);
var_eta = var(3,3);

Z1 = -1.64;
Z2 = 1.64;



q_ub = sol(1) + Z2*(var_q)^(0.5);
q_lb = sol(1) + Z1*(var_q)^(0.5);
beta_ub = sol(2) + Z2*(var_beta)^(0.5);
beta_lb = sol(2) + Z1*(var_beta)^(0.5);
eta_ub = sol(3) + Z2*(var_eta)^(0.5);
eta_lb = sol(3) + Z1*(var_eta)^(0.5);

ACI(1,1) = q_lb;
ACI(1,2) = q_ub;
ACI(1,3) = beta_lb;
ACI(1,4) = beta_ub;
ACI(1,5) = eta_lb;
ACI(1,6) = eta_ub;
ACI(1,7) = q_ub-q_lb;
ACI(1,8) = beta_ub-beta_lb;
ACI(1,9) = eta_ub-eta_lb;
    
    
 

