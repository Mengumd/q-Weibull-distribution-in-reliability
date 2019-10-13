global lifes;
global Cens;

theta = 0.2;
n = 100;
d = 4;
eta = 5*ones(1,d);
beta = ones(1,d);

rep = 10;
sols = zeros(rep,3);
objs = zeros(rep,1);

Cens = zeros(1,n);
lifes = zeros(1,n);

for ii = 1:n

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

    lifes(ii) = max(t);

end


[Obj, sol, output]=ABC_direct_throw_local_search_adaptive_qfrechet(1);


