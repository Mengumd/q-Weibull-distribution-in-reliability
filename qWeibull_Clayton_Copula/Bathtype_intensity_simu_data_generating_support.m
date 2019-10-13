global lifes;
global Cens;

theta = -0.2857;
d = 4;
eta = 5*ones(1,d);
beta = [2 2 2 2];

rep = 1;
sols = zeros(rep,3);
objs = zeros(rep,1);

for m = 1:rep
    
    Cens = zeros(1,100);
    lifes = zeros(1,100);

    t = zeros(1,d);
    R = ones(1,d);

    inslog = ( sum(R.^(-theta))  - d + 1)* (rand^(-theta)) + d - 1 -  (sum(R.^(-theta)) - R(1)^(-theta));
    t(1) = eta(1)*((1/theta)*log( inslog ) )^(1/beta(1));

    for gg = 2:d
        R(1:gg-1) = exp(-(t(1:gg-1)./eta(1:gg-1)).^beta(1:gg-1));
        R(gg:d) = 1;
        inslog = ( sum(R.^(-theta)) - d + 1)* (rand^(-theta/(1+(gg-1)*theta))) + d - 1 - (sum(R.^(-theta)) - R(gg)^(-theta));
        t(gg) = eta(gg)*((1/theta)*log(inslog))^(1/beta(gg));
    end

    lifes(1) = min(t);
    t_max = lifes(1);
    R = exp(-(t_max./eta).^beta);
    ii = 1;
    while(sum(R.^(-theta)) - d + 1 > 0.01)
        ii = ii + 1;                        
        t = zeros(1,d);
        R = exp(-(lifes(ii-1)./eta).^beta);
        inslog = ( sum(R.^(-theta))  - d + 1)* (rand^(-theta)) + d - 1 -  (sum(R.^(-theta))-R(1)^(-theta) );
        t(1) = eta(1)*((1/theta)*log(inslog ))^(1/beta(1));

        for jj = 2:d
            R(1:jj-1) = exp(-(t(1:jj-1)./eta(1:jj-1)).^beta(1:jj-1));
            R(jj:d) = exp(-(lifes(ii-1)./eta(jj:d)).^beta(jj:d));
            inslog = ( sum(R.^(-theta)) - d + 1)* (rand^(-theta/(1+(jj-1)*theta))) + d - 1 - (sum(R.^(-theta)) - R(jj)^(-theta));
            t(jj) = eta(jj)*((1/theta)*log(inslog ))^(1/beta(jj));
        end

        lifes(ii) = min(t);
        t_max = lifes(ii);
        R = exp(-(t_max./eta).^beta);
    end
                    
    n = ii;    
    Cens = Cens(1:n);
    lifes = lifes(1:n);

    [Obj, sol, output]=ABC_direct_throw_local_search_adaptive(1);
%     [Obj, sol, output]=ABC_direct_throw_local_search_adaptive_copula_theta_hr(1);

    sols(m,:) = sol;
    objs(m,:) = Obj;
    
    fprintf('sample number=%d \n',m);

end



y = sort(sols);
z = sort(objs);
% sol_lb = y(51,:);
% sol_ub = y(950,:);
% obj_lb = z(51);
% obj_ub = z(950);
% sol_lb = y(5,:);
% sol_ub = y(95,:);
% obj_lb = z(5);
% obj_ub = z(95);
    
save data_ApplicationExample_1_BCI_P_sample100
