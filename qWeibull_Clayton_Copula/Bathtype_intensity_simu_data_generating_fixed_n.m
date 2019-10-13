global lifes;
global Cens;

theta = -0.3333;
n = 100;
d = 4;
eta = 5*ones(1,d);
beta = [2 2 0.2 0.2];

rep = 10;
sols = zeros(rep,3);
objs = zeros(rep,1);

for m = 1:rep

    Cens = zeros(1,n);
    lifes = zeros(1,n);

    t = zeros(1,d);
    R = ones(1,d);

    inslog = ( sum(R.^(-theta))  - d + 1)* (rand^(-theta)) + d - 1 -  (sum(R.^(-theta)) - R(1)^(-theta));
    while(inslog > 1)
        inslog = ( sum(R.^(-theta))  - d + 1)* (rand^(-theta)) + d - 1 -  (sum(R.^(-theta)) - R(1)^(-theta));
    end
    t(1) = eta(1)*((1/theta)*log( inslog ) )^(1/beta(1));

    for gg = 2:d
        R(1:gg-1) = exp(-(t(1:gg-1)./eta(1:gg-1)).^beta(1:gg-1));
        R(gg:d) = 1;
        inslog = ( sum(R.^(-theta)) - d + 1)* (rand^(-theta/(1+(gg-1)*theta))) + d - 1 - (sum(R.^(-theta)) - R(gg)^(-theta));
        while(inslog >1)
            inslog = ( sum(R.^(-theta)) - d + 1)* (rand^(-theta/(1+(gg-1)*theta))) + d - 1 - (sum(R.^(-theta)) - R(gg)^(-theta));
        end
        t(gg) = eta(gg)*((1/theta)*log(inslog))^(1/beta(gg));
    end

    lifes(1) = min(t);

    for ii = 2:n
        t = zeros(1,d);
        R = exp(-(lifes(ii-1)./eta).^beta);
        inslog = ( sum(R.^(-theta))  - d + 1)* (rand^(-theta)) + d - 1 -  (sum(R.^(-theta))-R(1)^(-theta) );
        while(inslog > 1)
            inslog = ( sum(R.^(-theta))  - d + 1)* (rand^(-theta)) + d - 1 -  (sum(R.^(-theta))-R(1)^(-theta) );
        end
        t(1) = eta(1)*((1/theta)*log(inslog ))^(1/beta(1));

        for jj = 2:d
            R(1:jj-1) = exp(-(t(1:jj-1)./eta(1:jj-1)).^beta(1:jj-1));
            R(jj:d) = exp(-(lifes(ii-1)./eta(jj:d)).^beta(jj:d));
            inslog = ( sum(R.^(-theta)) - d + 1)* (rand^(-theta/(1+(jj-1)*theta))) + d - 1 - (sum(R.^(-theta)) - R(jj)^(-theta));
            while(inslog > 1)
                inslog = ( sum(R.^(-theta)) - d + 1)* (rand^(-theta/(1+(jj-1)*theta))) + d - 1 - (sum(R.^(-theta)) - R(jj)^(-theta));
            end
            t(jj) = eta(jj)*((1/theta)*log(inslog ))^(1/beta(jj));
        end

        lifes(ii) = min(t);
    end

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
