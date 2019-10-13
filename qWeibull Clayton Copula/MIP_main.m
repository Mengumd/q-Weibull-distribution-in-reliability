
global lifes;

load('LifeDataFile');

objMean = zeros(20,1);
objStd = zeros(20,1);
solMean = zeros(20,3);
solStd = zeros(20,3);
for ii = 1:20
    q = ES(ii,1);
    beta = ES(ii,2);
    eta = ES(ii,3);
    n = ES(ii,4);
    lifes = LifeData{ii};
    tic;
    [Obj, sol]=ABC_direct_throw_local_search_ahabc(1);
    toc;
    objMean(ii) = mean(Obj);
    objStd(ii) = std(Obj);
    solMean(ii,:) = sol;
    solStd(ii,:) = std(sol);
end
save data_MIP_qWeibull