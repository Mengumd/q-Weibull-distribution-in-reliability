
global lifes;

theta = -0.2;
d = 4;
eta = 5*ones(1,d);
beta = 2*ones(1,d);
n = 100000;        % sample number

ts = zeros(n,d);
dts = zeros(n,d);
ts_independent = zeros(n,d);
dts_independent = zeros(n,d);

lifes = zeros(1,n);

for m = 1:n

    t = zeros(1,d);
    F = ones(1,d);

    inslog = ( sum(F.^(-theta))  - d + 1)* (rand^(-theta)) + d - 1 -  (sum(F.^(-theta)) - F(1)^(-theta));
    t(1) = eta(1)*((1/theta)*log( inslog ) )^(-1/beta(1));

    for gg = 2:d
        F(1:gg-1) = exp(-(t(1:gg-1)./eta(1:gg-1)).^(-beta(1:gg-1)));
        F(gg:d) = 1;
        inslog = ( sum(F.^(-theta)) - d + 1)* (rand^(-theta/(1+(gg-1)*theta))) + d - 1 - (sum(F.^(-theta)) - F(gg)^(-theta));
        t(gg) = eta(gg)*((1/theta)*log(inslog))^(-1/beta(gg));
    end
    t = sort(t);
    ts(m,:) = t;
    dt = t - [0 t(1:end-1)];
    dts(m,:) = dt;

end

for m = 1:n

    t = zeros(1,d);
    for gg = 1:d
        u = rand;
        t(gg) = eta(gg)*(-log(u))^(-1/beta(gg));
    end
    t = sort(t);
    ts_independent(m,:) = t;
    dt = t - [0 t(1:end-1)];
    dts_independent(m,:) = dt;
end

tt = 0:0.5:50;
load_ratio = zeros(length(tt)-1,d);
for gg = 1:d
    for kk = 1:length(tt)-1
        dt_selected = dts(ts(:,gg)>tt(kk) & ts(:,gg)<tt(kk+1));        
        dt_independent_selected = dts_independent(ts_independent(:,gg)>tt(kk) & ts_independent(:,gg)<tt(kk+1));
        mean(dt_selected)./mean(dt_independent_selected);
        load_ratio(kk,gg) = mean(dt_selected)./mean(dt_independent_selected);
    end
end
figure;plot(load_ratio');

dt_mean = mean(dts);
dt_independent_mean = mean(dts_independent);
figure;plot(dt_mean);
figure;plot(dt_independent_mean);
figure;plot(dt_mean./dt_independent_mean);
%     lifes = lifes(1:n);
% 
%     [Obj, sol, output]=ABC_direct_throw_local_search_adaptive_qfrechet(1);
%     sols(r,:)=sol;
%     objs(r,:)=Obj;

save Data_qfrechet_Interval_simulation
