% load Data_intensity_ks_test

global lifes;
global Cens;

% Application Example
[TBF,~,~] = xlsread('bathtype_intensity.xlsx','Sheet7','A1:A44');
[Cens,~,~] = xlsread('bathtype_intensity.xlsx','Sheet7','B1:B44');

lifes = zeros(1,length(TBF));
lifes(1) = TBF(1);
for i = 2:length(TBF)
    lifes(i) = lifes(i-1) + TBF(i) ;
end

% [Obj, sol, output]=ABC_direct_throw_local_search_adaptive_SPLP(1);
% Application Example

t = lifes;
d = 6;
beta_un = [1.1823 0.5896 7.3851 2.213 1.3512 33.0909];
eta_un = [104.0686 9.2704 1395.8181 6638.3318 2810258.5216 72870.6967];
theta_un = 16.8148;
Copula_unknown = zeros(1,length(t));
h_unknown = zeros(1,length(t));
for ii = 1:length(t)
    tt = t(ii);
    R = exp(-(tt./eta_un).^beta_un);
    h = beta_un.*tt.^(beta_un-1)./eta_un.^beta_un;
    if sum(R.^(-theta_un)) == inf        % numerical stability
        Copula_unknown(ii) = min(R);
    else
        Copula_unknown(ii) = (sum(R.^(-theta_un)) - d + 1)^(-(1/theta_un));
    end
    
    h_unknown(ii) = sum((R.^(-theta_un)).*h)/(Copula_unknown(ii)^(-theta_un));
end


rep = 1000;

sol = [theta_un, beta_un, eta_un];
H_unknown = cell(1,rep);
H_unknown{1} = -log(Copula_unknown);

sols = zeros(rep,13);
objs = zeros(rep,1);
sols(1,:) = sol;

LifeData = cell(1,rep);
CensData = cell(1,rep);
LifeData{1} = lifes;
n = length(lifes);

D = zeros(rep,1);
Fn = zeros(1,n+1);
for j = 0:n
    Fn(j+1,1) = j / n ;
end

Num = H_unknown{1};
Num_Ratio = Num/Num(end);

Dmais = max(abs(Fn(2:n+1) - Num_Ratio));
Dmenos = max(abs(Fn(1:n) - Num_Ratio));
D(1) = max(Dmais, Dmenos);

for m = 2:rep
    
    Cens = zeros(1,n);
    lifes = zeros(1,n);
    
    beta = [1.1823 0.5896 7.3851 2.213 1.3512 33.0909];
    eta = [104.0686 9.2704 1395.8181 6638.3318 2810258.5216 72870.6967];
    theta = 16.8148;
    
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
    [tmin,comp_index] = min(t);
    Cens(1) = comp_index;    
    
    for ii = 2:n                    
        t = zeros(1,d);
        R = exp(-(lifes(ii-1)./eta).^beta);
        u = rand;
%         inslog = ( sum(R.^(-theta))  - d + 1)* (u^(-theta)) + d - 1 -  (sum(R.^(-theta))-R(1)^(-theta) );
        inslog1 = (sum(R.^(-theta)))* (u^(-theta) - 1);
        inslog2 = (- d + 1)* (u^(-theta)) + d - 1 + R(1)^(-theta);
        if inslog1 == inf %&& inslog2 ~= inf
                loginslog1 = -theta*log(min(R))+log((u^(-theta) - 1));
                t(1) = eta(1)*((1/theta)*loginslog1)^(1/beta(1));
            else
                %if inslog1 ~= inf && inslog2 ~= inf
                inslog = inslog1 + inslog2;
                t(1) = eta(1)*((1/theta)*log(inslog ))^(1/beta(1));
                %else 
                %    bug = 1;
                %    display('bug found');
                %end
        end
        
        for jj = 2:d
            R(1:jj-1) = exp(-(t(1:jj-1)./eta(1:jj-1)).^beta(1:jj-1));
            R(jj:d) = exp(-(lifes(ii-1)./eta(jj:d)).^beta(jj:d));
            
            u = rand;
            %inslog = ( sum(R.^(-theta)) - d + 1)* (u^(-theta/(1+(jj-1)*theta))) + d - 1 - (sum(R.^(-theta)) - R(jj)^(-theta));
            inslog1 = sum(R.^(-theta))* (u^(-theta/(1+(jj-1)*theta))-1);
            inslog2 = (- d + 1)* (u^(-theta/(1+(jj-1)*theta))) + d - 1 + R(jj)^(-theta);
            if inslog1 == inf %&& inslog2 ~= inf
                loginslog1 = -theta*log(min(R))+log((u^(-theta/(1+(jj-1)*theta))-1));
                t(jj) = eta(jj)*((1/theta)*loginslog1)^(1/beta(jj));
            else
                %if inslog1 ~= inf && inslog2 ~= inf
                inslog = inslog1 + inslog2;
                t(jj) = eta(jj)*((1/theta)*log(inslog ))^(1/beta(jj));
%                 else 
%                     bug = 1;
%                     display('bug found');
%                 end
            end
        end

        lifes(ii) = min(t);
        t_max = lifes(ii);
        R = exp(-(t_max./eta).^beta);
        [tmin,comp_index] = min(t);
        Cens(ii) = comp_index;
        
        
        if t_max == inf
            error = 1;
        end
        
        
        if sum(R.^(-theta)) - d + 1 < 0.001
            break;
        end
        
    end
    
    LifeData{m} = lifes;
    CensData{m} = Cens;
end

for m = 2:rep
    
    fprintf('sample number=%d \n',m);
    
    lifes = LifeData{m};
    Cens = zeros(1,length(lifes));
    
%     [Obj, sol, output]=ABC_direct_throw_local_search_adaptive(1);
%     [Obj, sol, output]=ABC_direct_throw_local_search_adaptive_copula_theta(1);
    [Obj, sol, output]=ABC_direct_throw_local_search_adaptive_copula_theta_hr(1);

    sols(m,:) = sol;
    objs(m,:) = Obj;
%     
%     [sol,obj]=fminsearch(@(x)objfun(x),[q,beta,eta]);
end

for m = 2:rep
    
    lifes = LifeData{m};
    t = lifes;
    sol = sols(m,:);
    
    theta = -sol(1);
    beta = sol(2:7);
    eta = sol(8:13);
    
    Copula_unknown = zeros(1,length(t));
    h_unknown = zeros(1,length(t));
    for ii = 1:length(t)
        tt = t(ii);
        R = exp(-(tt./eta).^beta);
        h = beta.*tt.^(beta-1)./eta.^beta;
        
        if sum(R.^(-theta)) == inf        % numerical stability
            Copula_unknown(ii) = min(R);
        else
            Copula_unknown(ii) = (sum(R.^(-theta)) - d + 1)^(-(1/theta));
        end
        h_unknown(ii) = sum((R.^(-theta)).*h)/(Copula_unknown(ii)^(-theta));
    end

    H_unknown{m} = -log(Copula_unknown);
    
    Fn = zeros(1,length(lifes)+1);
    for j = 0:length(lifes)
        Fn(j+1,1) = j / length(lifes) ;
    end
    
    Num = H_unknown{m};
    Num_Ratio = Num/Num(end);

    Dmais = max(abs(Fn(2:length(lifes)+1) - Num_Ratio));
    Dmenos = max(abs(Fn(1:length(lifes)) - Num_Ratio));
    D(m) = max(Dmais, Dmenos);
    
    fprintf('sample number=%d \n',m);

end

a = find(D > D(1)) ;
p = (size(a,1) + 1) / rep

figure;
t_d = LifeData{1};

Fn = zeros(1,length(t_d)+1);
for j = 0:length(t_d)
    Fn(j+1) = j;
end

stairs([0 t_d],Fn,'b-');
hold on;
t = 0:0.01:t_d(end);
R = (1-(1-q)*(t./eta).^beta).^((2-q)/(1-q));
H = -log(R);
plot (t, H,'k--');



figure;
stairs([0 lifes],Fn,'k-');
hold;
plot(t,H_known,'m-');
legend('Empirical','Copula Component Known');
title('Expected Number of Failures');
xlabel('time');
ylabel('H(t)');

figure;
plot(t,h_known,'m-');
legend('Copula Component Known');
title('Intensity Function');
xlabel('time');
ylabel('h(t)');

