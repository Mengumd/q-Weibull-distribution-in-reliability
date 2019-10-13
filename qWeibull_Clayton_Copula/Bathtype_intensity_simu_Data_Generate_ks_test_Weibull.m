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

% [Obj, sol, output]=ABC_direct_throw_local_search_adaptive_weibull(1);
% Application Example

beta = 0.9257;
eta = 40.8725;
sol = [beta,eta];

rep = 1000;
sols = zeros(rep,2);
objs = zeros(rep,1);
sols(1,:) = sol;

LifeData = cell(1,rep);
LifeData{1} = lifes;

R_weibull = cell(1,rep);
H_weibull = cell(1,rep);

n = length(lifes);
t = lifes;
R_weibull{1} = exp(-(t/eta).^beta);
H_weibull{1} = -log(R_weibull{1});

D = zeros(rep,1);
Fn = zeros(1,n+1);
for j = 0:n
    Fn(j+1,1) = j / n ;
end

Num = H_weibull{1};
Num_Ratio = Num/Num(end);

Dmais = max(abs(Fn(2:n+1) - Num_Ratio));
Dmenos = max(abs(Fn(1:n) - Num_Ratio));
D(1) = max(Dmais, Dmenos);

for m = 2:rep
    
    lifes = zeros(1,n);
    lifes(1) = qwblrnd(1,eta,beta,1);        
    for ii = 2:n
        cdt = exp(-(lifes(ii-1)/eta)^beta);
        t_rnd = qwblrnd_cdt(1,eta,beta,cdt);        
        lifes(ii) = t_rnd;
        
%         if(t_rnd>2317)
%             break;
%         end
        
    end
%     lifes(ii+1:end)=[];
    
    LifeData{m} = lifes;
    
    Cens = zeros(1,n);
    [obj, sol, output]=ABC_direct_throw_local_search_adaptive_weibull(1);
%     [sol,obj]=fminsearch(@(x)objfun(x),[q,beta,eta]);

    sols(m,:) = sol;
    objs(m,:) = obj;
    
    t = lifes;
    R_weibull{m} = exp(-(t/sol(2)).^sol(1));
    H_weibull{m} = -log(R_weibull{m});
    
    
    Fn = zeros(1,length(lifes)+1);
    for j = 0:length(lifes)
        Fn(j+1,1) = j / length(lifes) ;
    end
    
    Num = H_weibull{m};
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

legend('NHPP qWeibull','Empirical');
ylabel('H(t)');
xlabel('Time');
title('Cumulative Intensity Function');
axis([0 2000 0 50]);

















































figure;
stairs([0 lifes],Fn,'k-');
hold;
plot(t,H_qweibull,'b-');
legend('Empirical','qWeibull');
title('Expected Number of Failures');
xlabel('time');
ylabel('H(t)');


save Data_intensity_ks_test_1000
