function ApplicationExample_2_bootstrap_P_KS_test_weibull()
global lifes;
lifes = [0.058 0.070 0.090 0.105 0.113 0.121 0.153 0.159 0.224 0.421 0.570 0.596 0.618 0.834 1.019 1.104 1.497 2.027 2.234 2.372 2.433 2.505 2.690 2.877 2.879 3.166 3.455 3.551 4.378 4.872 5.085 5.272 5.341 8.952 9.188 11.399];

%% solve the mle of original data

Func_Paras.D = 2; %/*The number of parameters of the problem to be optimized*/
Func_Paras.ub = [100 100];%[1.9 10 max_eta]; %/*lower bounds of the parameters. */
Func_Paras.lb = [0 0];%[-10 0.00001 0.0001];%/*upper bound of the parameters.*/
Func_Paras.getNewFeasibleSolution = @()getNewFeasibleSolution(Func_Paras.lb,Func_Paras.ub,Func_Paras.D);
Func_Paras.anyObjFun = @negative_likelihood_func;
Func_Paras.isFeasible = @isFeasible;

ABC_Paras.NP=100; %/* The number of colony size (employed bees+onlooker bees)*/
ABC_Paras.FoodNumber=ABC_Paras.NP/2; %/*The number of food sources equals the half of the colony size*/
ABC_Paras.limit=150; %/*A food source which could not be improved through "limit" trials is abandoned by its employed bee*/
ABC_Paras.maxCycle=50000; %/*The number of cycles for foraging {a stopping criteria}*/
ABC_Paras.maxBestTrial = 1000;
ABC_Paras.maxnFcn = 2e5;
ABC_Paras.deltaObj = 1e-16;
ABC_Paras.isLocalAdded = 1;
ABC_Paras.method = 'AHABC';
ABC_Paras.C = 1;
ABC_Paras.Ns = 5;
ABC_Paras.runtime = 10;%/*Algorithm can be run many times in order to see its robustness*/

[obj, sol, ~]=ABC_ahabc_allFunc(Func_Paras,ABC_Paras);
solMean = mean(sol);
solstd = std(sol);
objMean = mean(obj);
objstd = std(obj);


save data_ApplicationMIP

%% define the MLE problems
function neg_loglik=negative_likelihood_func(paras)
global nFcn;
nFcn = nFcn+1;

global lifes;
t = lifes;

t0 =0;

beta = paras(1);
eta = paras(2);
q = 1;
neg_loglik_vector = -(log((2-q)*beta/(eta-t0)) + (beta-1)*log((t-t0)/(eta-t0)) + log(exp_q(-((t-t0)./(eta-t0)).^beta, q)));
neg_loglik = sum(neg_loglik_vector);


function b = isFeasible(paras)
b = paras(1)>0 && paras(2)>0;

function sol = getNewFeasibleSolution(lb,ub,D)
sol = rand(1,D).*(ub-lb) + lb;
while ~isFeasible(sol)
    sol = rand(1,D).*(ub-lb) + lb;
end

