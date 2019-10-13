% load Data_intensity_ks_test

global lifes;
global Cens;

% Application Example
[TBF,~,~] = xlsread('bathtype_intensity.xlsx','Sheet7','A1:A44');
[Cens,~,~] = xlsread('bathtype_intensity.xlsx','Sheet7','B1:AB44');

lifes = zeros(1,length(TBF));
lifes(1) = TBF(1);
for i = 2:length(TBF)
    lifes(i) = lifes(i-1) + TBF(i) ;
end

[Obj, sol, output]=ABC_direct_throw_local_search_adaptive(1);
% Application Example

q = sol(1);
beta = sol(2);
eta = sol(3);


% q = 0.9495;
% beta = 0.5652;
% eta = 15.1495; 
sol = [q,beta,eta];

rep = 1000;
sols = zeros(rep,3);
objs = zeros(rep,1);
sols(1,:) = sol;

LifeData = cell(1,rep);
LifeData{1} = lifes;

R_qweibull = cell(1,rep);
H_qweibull = cell(1,rep);

n = length(lifes);
t = lifes;
R_qweibull{1} = (1-(1-q)*(t./eta).^beta).^((2-q)/(1-q));
H_qweibull{1} = -log(R_qweibull{1});

D = zeros(rep,1);
Fn = zeros(1,n+1);
for j = 0:n
    Fn(j+1,1) = j / n ;
end

Num = H_qweibull{1};
Num_Ratio = Num/Num(end);

Dmais = max(abs(Fn(2:n+1) - Num_Ratio));
Dmenos = max(abs(Fn(1:n) - Num_Ratio));
D(1) = max(Dmais, Dmenos);

for m = 2:rep
    
    lifes = zeros(1,n);
    lifes(1) = qwblrnd(q,eta,beta,1);        
    for ii = 2:n
        cdt = (1-((1-q)*(lifes(ii-1)/eta)^beta))^((2-q)/(1-q));
        t_rnd = qwblrnd_cdt(q,eta,beta,cdt);        
        lifes(ii) = t_rnd;
        
%         if(t_rnd>2317)
%             break;
%         end
        
    end
%     lifes(ii+1:end)=[];
    
    LifeData{m} = lifes;
    
    Cens = zeros(1,n);
    [obj, sol, output]=ABC_direct_throw_local_search_adaptive(1);
%     [sol,obj]=fminsearch(@(x)objfun(x),[q,beta,eta]);

    sols(m,:) = sol;
    objs(m,:) = obj;
    
    t = lifes;
    R_qweibull{m} = (1-(1-sol(1))*(t./sol(3)).^sol(2)).^((2-sol(1))/(1-sol(1)));
    H_qweibull{m} = -log(R_qweibull{m});
    
    
    Fn = zeros(1,length(lifes)+1);
    for j = 0:length(lifes)
        Fn(j+1,1) = j / length(lifes) ;
    end
    
    Num = H_qweibull{m};
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
    Fn(j+1) = j/n;
end

stairs([0 t_d],Fn,'b-');
hold on;
t = 0:0.01:t_d(end);
R = (1-(1-q)*(t./eta).^beta).^((2-q)/(1-q));
H = -log(R);
H = H/H(end);
plot (t, H,'k--');

legend('Empirical Model','NHPP qWeibull Model');
ylabel('$\displaystyle\frac{H(t)}{H(t_{n})}$','interpreter','latex');
xlabel('Time');
title('Normalized Cumulative Intensity Function');
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
