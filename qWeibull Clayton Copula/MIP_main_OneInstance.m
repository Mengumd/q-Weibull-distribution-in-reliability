
global lifes;

LifeData=cell(1,50);
ES=[0.5 0.5 5 20;
    0.5 0.5 5 100;
    0.5 0.5 5 500;
    0.5 0.5 5 1000;
    1.5 0.5 5 50;
    1.5 0.5 5 100;
    1.5 0.5 5 500;
    1.5 0.5 5 1000;
    1 1 5 20;
    1 1 5 100;
    1 1 5 500;
    1 1 5 1000;
    0.5 1.5 5 20;
    0.5 1.5 5 100;
    0.5 1.5 5 500;
    0.5 1.5 5 1000;
    1.5 1.5 5 50;
    1.5 1.5 5 100;
    1.5 1.5 5 500;
    1.5 1.5 5 1000;];

N = 30;

for ii = 8  
    q = ES(ii,1);
    beta = ES(ii,2);
    eta = ES(ii,3);
    n = ES(ii,4);
    
    obj_qWeibull = zeros(N,1);
    sol_qWeibull = zeros(N,3);
    
    obj_Weibull = zeros(N,1);
    sol_Weibull = zeros(N,2);
    
    for jj = 1:N
        LifeData{jj} = qwblrnd(q,eta,beta,n);
        
        lifes = LifeData{jj};
        tic;
        [obj, sol]=ABC_direct_throw_local_search_ahabc(1);
        toc;

        obj_qWeibull(jj) = obj;
        sol_qWeibull(jj,:) = sol;
        
        tic;
        [obj, sol]=MIP_Weibull_OneInstance();
        toc;
        
        obj_Weibull(jj) = obj;
        sol_Weibull(jj,:) = sol;
        
    end
     
        
end

save data_MIP_qWeibull_OneInstance_SampleSize100_8