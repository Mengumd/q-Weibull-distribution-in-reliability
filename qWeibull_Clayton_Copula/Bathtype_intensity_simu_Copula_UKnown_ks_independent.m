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
beta = [0.6082 0.6082 156.2540 0.6082 2.6416 0.6082];
eta = [161.5754 54.3404 2302.1496 75.0654 788.0360 1998.8373];

Copula_unknown = zeros(1,length(t));
h_unknown = zeros(1,length(t));
for ii = 1:length(t)
    tt = t(ii);
    R = exp(-(tt./eta).^beta);
    Copula_unknown(ii) = prod(R);
end

rep = 1000;

sol = [beta, eta];
H_unknown = cell(1,rep);
H_unknown{1} = -log(Copula_unknown);

sols = zeros(rep,12);
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
    
    beta = [0.6082 0.6082 156.2540 0.6082 2.6416 0.6082];
    eta = [161.5754 54.3404 2302.1496 75.0654 788.0360 1998.8373];
    
    t = zeros(1,d);
    for gg = 1:d
        u = rand;
        t(gg) = eta(gg)*((-log(u)) )^(1/beta(gg));
    end
    lifes(1) = min(t);
    [tmin,comp_index] = min(t);
    Cens(1) = comp_index;    
    
    for ii = 2:n
        t = zeros(1,d);
        for jj = 1:d
            cdt = exp(-(lifes(ii-1)/eta(jj)).^beta(jj));
            u = rand*cdt;
            t(jj) = eta(jj)*((-log(u)) )^(1/beta(jj));
        end
        
        lifes(ii) = min(t);
        [tmin,comp_index] = min(t);
        Cens(ii) = comp_index;
        
%         if(t_rnd>2317)
%             break;
%         end
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
    [Obj, sol, output]=copula_theta_unknown_independent(1);

    sols(m,:) = sol;
    objs(m,:) = Obj;
%     
%     [sol,obj]=fminsearch(@(x)objfun(x),[q,beta,eta]);
end

for m = 2:rep
    
    lifes = LifeData{m};
    t = lifes;
    sol = sols(m,:);
    
    beta = sol(1:6);
    eta = sol(7:12);
    
    Copula_unknown = zeros(1,length(t));
    for ii = 1:length(t)
        tt = t(ii);
        R = exp(-(tt./eta).^beta);
        Copula_unknown(ii) = prod(R);
    end
    
    sol = [beta, eta];
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

