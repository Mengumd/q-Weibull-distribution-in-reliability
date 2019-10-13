alpha_1 = 11.89;
beta_1 = 0.603;
alpha_2 = 912;
beta_2 = 3.211;

global lifes;
t = lifes;

for ii = 1:length(t)
    if Cens(ii) == 0
        loglik_vector(ii) = log((beta_1/alpha_1)*(t(ii)/alpha_1).^(beta_1-1) + (beta_2/alpha_2)*(t(ii)/alpha_2).^(beta_2-1)) ;
    end
end

loglik = sum(loglik_vector) - ((t(end)/alpha_1)^beta_1 + (t(end)/alpha_2)^beta_2)

% ((lifes(end)/alpha_1)^beta_1) + (lifes(end)/alpha_2)^beta_2

% -((2-q)*log(1-(1-q)*(lifes(end)/eta)^beta))/(1-q)