[TBF,~,~] = xlsread('bathtype_intensity.xlsx','Sheet7','A1:A44');
[Cens,~,~] = xlsread('bathtype_intensity.xlsx','Sheet7','B1:B44');

lifes = zeros(1,length(TBF));
lifes(1) = TBF(1);
for i = 2:length(TBF)
    lifes(i) = lifes(i-1) + TBF(i) ;
end

t = lifes;
for ii = 1:length(t)
    if Cens(ii) == 0
        loglik_vector(ii) = log(0.603*t(ii).^(0.603-1)/11.89^0.603 + 3.211*t(ii).^(3.211-1)/912^3.211);
    end
end

loglik = sum(loglik_vector) - ((t(end)/11.89).^0.603 + (t(end)/912).^3.211)


