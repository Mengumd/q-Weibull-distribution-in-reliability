global lifes;
global Cens;

eta = 5;
beta = 2;
theta = 1;
n = 5;
d = 2;
a = 1/theta;

Cens = zeros(1,n);
lifes = zeros(1,n);

x = rand(1,d);
u = (1 - log(x)/gamrnd(a,1)).^(-1/theta);
t = eta*(-log(u)).^(1/beta);
lifes(1) = min(t);

for ii = 2:n
    x = rand(1,d);
    u = (1 - log(x)/gamrnd(a,1)).^(-1/theta);    
    t = eta*(-log(u)).^(1/beta);
    
    while min(t) <= lifes(ii-1)
        x = rand(1,d);
        u = (1 - log(x)/gamrnd(a,1)).^(-1/theta);    
        t = eta*(-log(u)).^(1/beta);
    end
    
    lifes(ii) = min(t);
end



[Obj, sol, output]=ABC_direct_throw_local_search_adaptive(1);


LifeData=cell(1,100);

q = sol(1);
beta = sol(1,2);
eta = sol(1,3);
n = 44;
for jj = 2:100
    LifeData{jj} = qwblrnd(q,eta,beta,n);
end

BCI_P_sol = zeros(1,9);
BCI_P_obj = zeros(1,3);

solutions = zeros(100,3);
objection = zeros(100,1);

solutions(1,:) = sol;

for jj = 2:100
    lifes = LifeData{jj};
    Cens = zeros(44,1);
    [Obj, sol, output]=ABC_direct_throw_local_search_adaptive(1);
    solutions(jj,:) = sol;
    objection(jj) = Obj;
         
    fprintf('sample number=%d \n',jj);
end
y = sort(solutions);
z = sort(objection);
% sol_lb = y(51,:);
% sol_ub = y(950,:);
% obj_lb = z(51);
% obj_ub = z(950);
sol_lb = y(5,:);
sol_ub = y(95,:);
obj_lb = z(5);
obj_ub = z(95);
    
BCI_P_sol(1) = sol_lb(1);
BCI_P_sol(2) = sol_ub(1);
BCI_P_sol(3) = sol_lb(2);
BCI_P_sol(4) = sol_ub(2);
BCI_P_sol(5) = sol_lb(3);
BCI_P_sol(6) = sol_ub(3);
BCI_P_sol(7) = sol_ub(1)-sol_lb(1);
BCI_P_sol(8) = sol_ub(2)-sol_lb(2);
BCI_P_sol(9) = sol_ub(3)-sol_lb(3);
BCI_P_obj(1) = obj_lb;
BCI_P_obj(2) = obj_ub;
BCI_P_obj(3) = obj_ub-obj_lb;

save data_ApplicationExample_1_BCI_P_sample100
