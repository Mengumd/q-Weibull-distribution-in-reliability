function [GlobalMins,solutions, plotOutPut]=ABC_ahabc_GRP(isLocalAdded)
%/* Control Parameters of ABC algorithm*/

NP=100; %/* The number of colony size (employed bees+onlooker bees)*/
FoodNumber=NP/2; %/*The number of food sources equals the half of the colony size*/
limit=150; %/*A food source which could not be improved through "limit" trials is abandoned by its employed bee*/
maxCycle=50000; %/*The number of cycles for foraging {a stopping criteria}*/
maxBestTrial = 1000;
maxnFcn = 5e5;
deltaObj = 1e-16;
Nscout = 0;

global nFcn;
global lifes;

D=4; %/*The number of parameters of the problem to be optimized*/
max_eta = mean(lifes);
ub=[1.9 10 max_eta 1]; %/*lower bounds of the parameters. */
lb=[-10 0.00001 0.0001 0.00001];%/*upper bound of the parameters.*/

runtime = 10;%/*Algorithm can be run many times in order to see its robustness*/

plotOutPut = cell(1,runtime);


GlobalMins=zeros(1,runtime);
solutions = [];

for r=1:runtime
  
    fprintf('run time=%d \n',r);
%% initialize
max_t = max(lifes);
nFcn = 0;
Nscout = 0;

Foods = zeros(FoodNumber,4);
  for i = 1:FoodNumber
      Foods(i,:) = getNewFeasibleSolution(lb, ub);
  end

ObjVal = zeros(1,FoodNumber);
    for iFoods=1:FoodNumber
        ObjVal(iFoods)=objfun(Foods(iFoods,:));
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

        ObjValSol=objfun(sol);
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
       
        ObjValSol=objfun(sol);
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
    
    ObjValSol=objfun(sol);
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
    [sol, ObjValSol, exitFlag]  = fminsearch(@(x)objfun((x)),Foods(indb,:),options);
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
% sol = fminsearch(@(x)objfun((x)),GlobalParams,options);
% if (isFeasible(sol))
%     ObjValSol=objfun(sol);
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

    b = 1;
    
    global lifes;
    t = lifes;
    
if ( sol(1) >= 2 )
    b = 0;
end

if (sol(2) <= 0)
    b = 0;
end

if (sol(3) <= 0)
    b = 0;
end

if (sol(4) < 0)
    b = 0;
end

if (sol(4) > 1)
    b = 0;
end

if ((1-(1-sol(1))*(t(1)/sol(3))^sol(2)) <= 0)
    b = 0;
end

for ii = 2:length(t)
    if((1-(1-sol(1))*((t(ii)+sol(4)*sum(t(1:ii-1)))/sol(3))^sol(2)) <= 0)
        b = 0;
    end
end

for ii = 2:length(t)
    if ((1-(1-sol(1))*(sol(4)*sum(t(1:ii-1))/sol(3))^sol(2)) <= 0)
        b = 0;
    end
end



function sol = getNewFeasibleSolution(lb, ub)
D = 4;
Foods = rand(1,D).* (ub-lb) + lb;
while(~isFeasible(Foods))
    Foods = rand(1,D).* (ub-lb) + lb;
end
sol = Foods;



