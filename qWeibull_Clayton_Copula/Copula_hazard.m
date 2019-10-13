d = 10;
theta = -0.5;
beta = 5*ones(1,d);
eta = 5*ones(1,d);
t = 0:0.01:10;
hr = zeros(1,length(t));

for ii = 1:length(t)
    R = exp(-(t(ii)./eta).^beta);
    h = (beta./eta).*(t(ii)./eta).^(beta-1);
    C = (sum(R.^theta) - d + 1)^(1/theta);
    hr(ii) = sum((R.^theta).*h)/(C^theta);
end

plot(t,hr);

t_max = eta./((theta/(1-theta)).^(1./beta));
t_max_copula = ((-(1/theta)*log(1-(1/d))).^(1./beta)).*eta;