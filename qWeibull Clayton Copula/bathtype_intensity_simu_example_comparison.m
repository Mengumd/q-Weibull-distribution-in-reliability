
[TBF,~,~] = xlsread('bathtype_intensity.xlsx','Sheet7','A1:A44');
[Cens,~,~] = xlsread('bathtype_intensity.xlsx','Sheet7','B1:B44');

lifes = zeros(1,length(TBF));
lifes(1) = TBF(1);
for i = 2:length(TBF)
    lifes(i) = lifes(i-1) + TBF(i) ;
end

Fn = zeros(1,length(lifes)+1);
for j = 0:length(lifes)
    Fn(j+1) = j;
end

t = 0:2317;
% t = 0:2200;
d = 6;
q = 0.9495;
beta_q = 0.5652;
eta_q = 15.1495;
R_qweibull = (1-(1-q)*(t./eta_q).^beta_q).^((2-q)/(1-q));
h_qweibull = ((2-q)*beta_q*t.^(beta_q-1)/eta_q^beta_q)./(1-(1-q)*(t/eta_q).^beta_q);

beta_w = 0.9257;
eta_w = 40.8725;
R_weibull = exp(-(t./eta_w).^beta_w);
h_weibull = beta_w*t.^(beta_w-1)/eta_w^beta_w;


beta = [0.7722 0.4395 3.6608 0.8612 1.1892 2.3211];
eta = [114.658 20.3106 1937.3287 770.2656 687.9873 1913.7023];
theta = -0.0245;
Copula_known = zeros(1,length(t));
h_known = zeros(1,length(t));
for ii = 1:length(t)
    tt = t(ii);
    R = exp(-(tt./eta).^beta);
    h = beta.*tt.^(beta-1)./eta.^beta;
    Copula_known(ii) = (sum(R.^(-theta)) - d + 1)^(-(1/theta));
    h_known(ii) = sum((R.^(-theta)).*h)/(Copula_known(ii)^(-theta));
end


beta_ki = [0.9071 0.5142 4.3294 1.0629 1.4752 2.8359];
eta_ki = [126.32 26.3187 1682.14 628.7710 619.5290 1572.83];
Copula_known_ind = zeros(1,length(t));
h_known_ind = zeros(1,length(t));
for ii = 1:length(t)
    tt = t(ii);
    R = exp(-(tt./eta_ki).^beta_ki);
    h = beta_ki.*tt.^(beta_ki-1)./eta_ki.^beta_ki;
    Copula_known_ind(ii) = prod(R);
    h_known_ind(ii) = sum(h);
end



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
        [rmin,index] = min(R);
        Copula_unknown(ii) = rmin;
        h_unknown(ii) = h(index);
    else
        Copula_unknown(ii) = (sum(R.^(-theta_un)) - d + 1)^(-(1/theta_un));
        h_unknown(ii) = sum((R.^(-theta_un)).*h)/(Copula_unknown(ii)^(-theta_un));
    end
    
end


beta_ui = [0.6082 0.6082 156.2540 0.6082 2.6416 0.6082];
eta_ui = [161.5754 54.3404 2302.1496 75.0654 788.0360 1998.8373];
Copula_unknown_ind = zeros(1,length(t));
h_unknown_ind = zeros(1,length(t));
for ii = 1:length(t)
    tt = t(ii);
    R = exp(-(tt./eta_ui).^beta_ui);
    h = (beta_ui./eta_ui).*(tt./eta_ui).^(beta_ui-1);
    Copula_unknown_ind(ii) = prod(R);
    h_unknown_ind(ii) = sum(h);
end


H_qweibull = -log(R_qweibull);
H_weibull = -log(R_weibull);
H_known = -log(Copula_known);
H_known_ind = -log(Copula_known_ind);
H_unknown = -log(Copula_unknown);
H_unknown_ind = -log(Copula_unknown_ind);
H_sps = (t/11.89).^0.603 + (t/912).^3.211;
h_sps = 0.603*t.^(0.603-1)/11.89^0.603 + 3.211*t.^(3.211-1)/912^3.211;


figure;
stairs([0 lifes],Fn,'k-','LineWidth',1.5);
hold;
plot(t,H_sps,'g-',t,H_qweibull,'b-',t,H_weibull,'r-',t,H_known,'m-',t,H_unknown,'k-.',t,H_known_ind,'c-',t,H_unknown_ind,'r-.','LineWidth',1.5);
legend('Empirical','S-PLP','qWeibull','Weibull','Copula Component Known','Copula Component Unknown','Independent Component Known','Independent Component Unknown');
title('Expected Number of Failures');
xlabel('time');
ylabel('H(t)');

figure;
plot(t,h_sps,'g-',t,h_qweibull,'b-',t,h_weibull,'r-',t,h_known,'m-',t,h_unknown,'k-.',t,h_known_ind,'c-',t,h_unknown_ind,'r-.','LineWidth',1);
legend('S-PLP','qWeibull','Weibull','Copula Component Known','Copula Component Unknown','Independent Component Known','Independent Component Unknown');
title('Intensity Function');
xlabel('time');
ylabel('h(t)');

axis([0 2500 0 0.05])






