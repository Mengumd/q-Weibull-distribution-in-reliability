t=0:0.01:10;
theta = 0.5; 
R=(max(1e-10,(exp(-1*t.^0.5*theta)+exp(-2*t.^0.3*theta)-1))).^(1/ theta);
hrf=-diff(R)./diff(t)./R(1:end-1);
figure;
plot(hrf)