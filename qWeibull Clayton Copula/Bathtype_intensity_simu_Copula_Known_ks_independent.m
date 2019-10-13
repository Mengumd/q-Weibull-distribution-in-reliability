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
beta = [0.9071 0.5142 4.3294 1.0629 1.4752 2.8359];
eta = [126.32 26.3187 1682.14 628.7710 619.5290 1572.83];

Copula_known = zeros(1,length(t));
for ii = 1:length(t)
    tt = t(ii);
    R = exp(-(tt./eta).^beta);
    Copula_known(ii) = prod(R);
end


rep = 1000;

sol = [beta, eta];
H_known = cell(1,rep);
H_known{1} = -log(Copula_known);

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

Num = H_known{1};
Num_Ratio = Num/Num(end);

Dmais = max(abs(Fn(2:n+1) - Num_Ratio));
Dmenos = max(abs(Fn(1:n) - Num_Ratio));
D(1) = max(Dmais, Dmenos);

for m = 2:rep
    
    Cens = zeros(1,n);
    lifes = zeros(1,n);
    
    beta = [0.7722 0.4395 3.6608 0.8612 1.1892 2.3211];
    eta = [114.658 20.3106 1937.3287 770.2656 687.9873 1913.7023];
    
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
%     lifes(ii+1:end)=[];

    LifeData{m} = lifes;
    CensData{m} = Cens;
end
  

for m = 2:rep
    
    fprintf('sample number=%d \n',m); 
    
    lifes = LifeData{m};
    Cens = CensData{m};
    
    [Obj, sol, output]=ABC_direct_throw_local_search_adaptive_copula_theta_independent(1);

    sols(m,:) = sol;
    objs(m,:) = Obj;
    
%     [sol,obj]=fminsearch(@(x)objfun(x),[q,beta,eta]);
end
    

for m = 2:rep
    
    lifes = LifeData{m};
    t = lifes;
    beta = sol(1:6);
    eta = sol(7:12);
    
    
    Copula_known = zeros(1,length(t));
    for ii = 1:length(t)
        tt = t(ii);
        R = exp(-(tt./eta).^beta);
        Copula_known(ii) = prod(R);
    end
    
    sol = [beta, eta];
    H_known{m} = -log(Copula_known);
    
    Fn = zeros(1,length(lifes)+1);
    for j = 0:length(lifes)
        Fn(j+1,1) = j / length(lifes) ;
    end
    
    Num = H_known{m};
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

