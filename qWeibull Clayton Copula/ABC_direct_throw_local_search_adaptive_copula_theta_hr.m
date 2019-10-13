function [GlobalMins,solutions, plotOutPut]=ABC_direct_throw_local_search_adaptive_copula_theta_hr(isLocalAdded)
%/* Control Parameters of ABC algorithm*/

NP=100; %/* The number of colony size (employed bees+onlooker bees)*/
FoodNumber=NP/2; %/*The number of food sources equals the half of the colony size*/
limit=150; %/*A food source which could not be improved through "limit" trials is abandoned by its employed bee*/
maxCycle=50000; %/*The number of cycles for foraging {a stopping criteria}*/
maxBestTrial = 1000;
maxnFcn = 7e5;
deltaObj = 1e-16;
Nscout = 0;

global nFcn;
global lifes;
% global Para_beta;
% global Para_eta;

D=13; %/*The number of parameters of the problem to be optimized*/
% ub = 1;
% lb = 0.001;
ub=[-10 5 5 5 5 5 5 100 100 100 100 100 100]; %/*upper bound of the parameters. */
lb=[-20 0 0 0 0 0 0 0 0 0 0 0 0];%/*lower bound of the parameters.*/
% ub=[1 3 3 3 3 3 3 100 100 2000 1000 1000 2000] ; %/*lower bounds of the parameters. */
% lb=[0.001 0 0 0 0 0 0 0 0 0 0 0 0];%/*upper bound of the parameters.*/
% ub=[2 10 10 10 10 10 10 100 100 100 100 100 100] ; %/*lower bounds of the parameters. */
% lb=[0.001 0 0 0 0 0 0 0 0 0 0 0 0];%/*upper bound of the parameters.*/
% ub=[1 0.6 20]; %/*lower bounds of the parameters. */
% lb=[0.9 0.5 10];%/*upper bound of the parameters.*/


runtime = 1;%/*Algorithm can be run many times in order to see its robustness*/

plotOutPut = cell(1,runtime);


GlobalMins=zeros(1,runtime);
solutions = [];

for r=1:runtime
  
    fprintf('run time=%d \n',r);
%% initialize
nFcn = 0;
Nscout = 0;

Foods = zeros(FoodNumber,D);
  for i = 1:FoodNumber
      Foods(i,:) = getNewFeasibleSolution(lb, ub);
  end

ObjVal = zeros(1,FoodNumber);
    for iFoods=1:FoodNumber
        ObjVal(iFoods)=objfun_copula_theta_hr(Foods(iFoods,:));
    end
Fitness=calculateFitness(ObjVal);

trial=zeros(1,FoodNumber);

%/*The best food source is memorized*/
BestInd=find(ObjVal==min(ObjVal));
BestInd=BestInd(end);
GlobalMin=ObjVal(BestInd);
GlobalParams=Foods(BestInd,:);

iter=1;
bestTrial = 0;
preGlobalMin = GlobalMin;
while (nFcn <= maxnFcn)%(iter <= maxCycle)),
%% %%%%%%% EMPLOYED BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:(FoodNumber)
        
        %/*The parameter to be changed is determined randomly*/
        Param2Change=fix(rand*D)+1;
        
        %/*A randomly chosen solution is used in producing a mutant solution of the solution i*/
        neighbour=fix(rand*(FoodNumber))+1;
       
        %/*Randomly selected solution must be different from the solution i*/        
            while(neighbour==i)
                neighbour=fix(rand*(FoodNumber))+1;
            end;
        
       sol=Foods(i,:);
       %  /*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
       sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2;
               
       if (~isFeasible(sol))
           continue;
       end

        ObjValSol=objfun_copula_theta_hr(sol);
        FitnessSol=calculateFitness(ObjValSol);
               
        if (FitnessSol > Fitness(i))
            Foods(i,:)=sol;
            Fitness(i)=FitnessSol;
            ObjVal(i)=ObjValSol;
            trial(i)=0;
        else
            trial(i)=trial(i)+1; %/*if the solution i can not be improved, increase its trial counter*/
       end;
                 
    end;
%%%%%%%%%%%%%%%%%%%%%%%% CalculateProbabilities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prob=(0.9.*abs(Fitness)./max(abs(Fitness)))+0.1;
  
%% %%%%%%%%%%%%%%%%%%%%%% ONLOOKER BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=1;
t=0;
while(t<FoodNumber)
    if(rand<prob(i))
        t=t+1;
        
        Param2Change=fix(rand*D)+1;
        neighbour=fix(rand*(FoodNumber))+1;
        
            while(neighbour==i)
                neighbour=fix(rand*(FoodNumber))+1;
            end;
        
       sol=Foods(i,:);
       sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2;
        
       if (~isFeasible(sol))
           continue;
       end
       
        ObjValSol=objfun_copula_theta_hr(sol);
        FitnessSol=calculateFitness(ObjValSol);
       
       if (FitnessSol > Fitness(i))
            Foods(i,:)=sol;
            Fitness(i)=FitnessSol;
            ObjVal(i)=ObjValSol;
            trial(i)=0;
        else
            trial(i)=trial(i)+1;
       end;
    end;
    
    i=i+1;
    if (i==(FoodNumber)+1) 
        i=1;
    end;   
    
end; 

%/*The best food source is memorized*/
         indb=find(ObjVal==min(ObjVal));
         indb=indb(end);
         if (ObjVal(indb)<GlobalMin)
             GlobalMin=ObjVal(indb);
             GlobalParams=Foods(indb,:);
         end
         
         if (preGlobalMin > GlobalMin)
             bestTrial = 0; % reset bestTrial
         else 
             bestTrial = bestTrial + 1;
         end;
         
% stop critirion
%          dObj = preGlobalMin - GlobalMin;
%          if (((dObj ~= 0) && (dObj < deltaObj)) || (bestTrial > maxBestTrial))
%              fprintf('best trial is : %d \n',bestTrial);
%              break;c
%          else
%              preGlobalMin = GlobalMin;
%          end
                  
%% %%%%%%%%%% SCOUT BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind=find(trial==max(trial));
ind=ind(end);
if (trial(ind)>limit)
    trial(ind)=0;
    
    sol = getNewFeasibleSolution(lb, ub);
    Nscout = Nscout + 1;
    
    ObjValSol=objfun_copula_theta_hr(sol);
    FitnessSol=calculateFitness(ObjValSol);
    Foods(ind,:)=sol;
    Fitness(ind)=FitnessSol;
    ObjVal(ind)=ObjValSol;
    fprintf('scout happens at food %d\n',ind);
end;

%% local search

if (isLocalAdded)
    FitnessGlobalMin = calculateFitness(GlobalMin);
    fprintf('Nscout = %d ', Nscout);
    C = 1;
    Ns = C*Nscout*limit;
    options = optimset('TolX',1e-9,'MaxIter',1e10,'MaxFunEvals',Ns,'Display','off');
    [sol, ObjValSol, exitFlag]  = fminsearch(@(x)objfun_copula_theta_hr((x)),Foods(indb,:),options);
    if (isFeasible(sol))
        FitnessSol=calculateFitness(ObjValSol);
        if (FitnessSol > Fitness(indb))
            ObjVal(indb) = ObjValSol;
            Foods(indb,:) = sol;
        end
    end
end
% if (exitFlag == 1)
%     GlobalMin = ObjValSol;
%     GlobalParams = sol;
%     fprintf('stop due to local search complete.\n');
%     break;
% end

% population  property
% fprintf('number of infeasible %d\n',sum(ObjVal == inf));

fprintf('iter=%d ObjVal=%g\n',iter,GlobalMin);
iter=iter+1;

plotOutPut{r}(end+1,:) = [nFcn, GlobalMin, mean(ObjVal), GlobalParams, Nscout]';
end % End of ABC

% %% final local search
% FitnessGlobalMin = calculateFitness(GlobalMin);
% 
% options = optimset('TolX',1e-5,'MaxIter',1e+10,'MaxFunEvals',1e+10);
% sol = fminsearch(@(x)objfun_copula_theta_hr((x)),GlobalParams,options);
% if (isFeasible(sol))
%     ObjValSol=objfun_copula_theta_hr(sol);
%     FitnessSol=calculateFitness(ObjValSol);
%     if (FitnessSol > FitnessGlobalMin)
%          GlobalMin=ObjValSol;
%          GlobalParams=sol;
%     end
% end

GlobalMins(r)=GlobalMin;
solutions(end+1,:)=GlobalParams;
end; %end of runs
save all

function fFitness=calculateFitness(fObjV)
fFitness=zeros(size(fObjV));
ind=find(fObjV>=0);
fFitness(ind)=1./(fObjV(ind)+1);
ind=find(fObjV<0);
fFitness(ind)=1+abs(fObjV(ind));

function b = isFeasible(sol)
if ( sol(1)<1 && sol(1)~=0 && min(sol(2:13))>0 )
% if ( sol(1)<1 && sol(1)~=0 )
    b = 1;
else
    b = 0;
end


function sol = getNewFeasibleSolution(lb, ub)
Foods = rand * (ub-lb) + lb;
sol = Foods;



