function [fitnessBest, popBest, fitnessCurrentBest_GA,popCurrentBest_GA] = ...
    ga(cost, sInitial,UB,LB, popSize, maxGeneration, selectionMethod, ...
    pCrossover, nElitism, numV, absTol, stallGeneration, imagineY)
%This RealGA is a minimization problem
%%Initialization of popultion and fitnessvalue
%Intial value of different patition

% xiaodong.zhang@u.nus.edu
% Last update August 05, 2021
% MATLAB version R2020b
%--------------------------------------------------------------------------------
popChildren = zeros(popSize-nElitism,numV);
fitnessChildren = zeros(1, popSize-nElitism);

popCurrent = sInitial;


popCurrent_GA = zeros(maxGeneration*popSize,numV);
fitnessCurrent_GA = zeros(1,maxGeneration*popSize);

popCurrentBest_GA = zeros(maxGeneration,numV);
fitnessCurrentBest_GA = zeros(1,maxGeneration);
fitnessCurrent = zeros(1,popSize);

for k = 1:popSize %Evaluate the fitness value of initial population %1
    fitnessCurrent(k) = cost(popCurrent(k,:));
end

%do maxGeneration time iteration, each iteration, the popCurrent = popChildren + popElitism

countG = 1; %the current number of generation
flag1 = 1; %for checking generation number
flag2 = 1; %for checking absoluate tol

while flag1 && flag2
    %%popElitism and fitnessChildren
    [fitnessCurrent_sorted, fitnessCurrent_Index] = sort(fitnessCurrent);
    fitnessElitism = fitnessCurrent_sorted(1:nElitism);
    popElitism = popCurrent(fitnessCurrent_Index(1:nElitism),:);
    
    for m = 1:popSize-nElitism  %%2main procedure, crossover, mutation and evaluation
        popCrossover = GAcrossOver(popCurrent, fitnessCurrent, selectionMethod, pCrossover);
        popMutation = GAmutation(popCrossover, numV, countG, maxGeneration);
        popChildren(m,:) = GAboundCheck(popMutation, UB, LB);
        fitnessChildren(m) =  cost(popChildren(m,:));
    end
    %%Return the current Population, current fitness, best Population and
    %%fitness among each generation
    
    %% population and fitness for current generation
    popCurrent = [popElitism;popChildren];% the populuation for current generation
    fitnessCurrent = [fitnessElitism, fitnessChildren];
    [fitnessCurrentBest, fCB_Index] = min(fitnessCurrent);
    popCurrentBest = popCurrent(fCB_Index,:);
    
    %% population and fitness for all the generations
    popCurrent_GA(1+(countG-1)*popSize:countG*popSize,:) = popCurrent;
    fitnessCurrent_GA(1+(countG-1)*popSize:countG*popSize) = fitnessCurrent;
    popCurrentBest_GA(countG,:) = popCurrentBest;
    fitnessCurrentBest_GA(countG) = fitnessCurrentBest;
    
    popBest = popCurrentBest_GA(countG,:);
    fitnessBest = fitnessCurrentBest_GA(countG);
    
    flag1 =  countG < maxGeneration;
    if countG <= stallGeneration || fitnessCurrentBest_GA(countG)==imagineY
        flag2 = 1;
    else
        if countG > 500 && fitnessCurrentBest_GA(countG)==imagineY
            flag2 = 0;
        else
            flag2 = abs(fitnessCurrentBest_GA(countG)-fitnessCurrentBest_GA(countG-stallGeneration)) > absTol;
        end
    end
    countG = countG+1;
end

end

function popCrossover = GAcrossOver(popCurrent, fitnessCurrent, selectionMethod, pCrossover)
% for each crossover, there is only one children

switch selectionMethod
    case 'Tournament'
        selection = @GAselectionT;
    case 'Roulette'
        selection = @GAselectionR;
end


popSelection1 = selection(popCurrent, fitnessCurrent);
popSelection2 = selection(popCurrent, fitnessCurrent);

if rand < pCrossover
    w = rand(1);  %randomly generate the weight
    popCrossover = w*popSelection1+(1-w)*popSelection2;
else
    if randi(2,1)==1
        popCrossover = popSelection1;
    else
        popCrossover = popSelection2;
    end
end

end


function popSelection = GAselectionR(popCurrent, fitnessCurrent)

imagineY=1000;

popCurrent = popCurrent(find(fitnessCurrent<imagineY),:);
fitnessCurrent = fitnessCurrent(find(fitnessCurrent<imagineY));
popSize = numel(fitnessCurrent);


[fitnessSorted, fitnessIndex] = sort(fitnessCurrent);
dis = abs(fitnessSorted-fitnessSorted(1)); %the distance to max fitness value
sum_dis = peak2peak(fitnessSorted);
pRoulette = dis./sum_dis;

randV = rand(1);

for k = 1:popSize-1 %the loop will only excute to popsize-1
    if pRoulette(k)<randV && randV<=pRoulette(k+1)
        popSelection = popCurrent(fitnessIndex(k),:);
        break;
    end
end



end


function popSelection = GAselectionT(popCurrent, fitnessCurrent)
%selection_T uses Tournament selection to selection one solution from
%current population

popSize = numel(fitnessCurrent);

rand1 = randi(popSize,1); %A random number between[1,popSize]
rand2 = randi(popSize,1);
pSelection1 = popCurrent(rand1,:);
pSelection2 = popCurrent(rand2,:);
fitness1 = fitnessCurrent(rand1);
fitness2 = fitnessCurrent(rand2);
if fitness1 <= fitness2
    popSelection = pSelection1;
else
    popSelection = pSelection2;
end

end


function popMutation = GAmutation(popCrossover, numV, numGeneration, maxGeneration)
% generate on mutation population from one crossover population

popMutation = zeros(1,numV);
pMutation = exp(-numGeneration/maxGeneration);

for m = 1:numV
    randv = rand;
    if randv < pMutation
        popMutation(m) = popCrossover(m)+ normrnd(0,abs(randv*popCrossover(m)));
    else
        popMutation(m) = popCrossover(m);
    end
end

end

function x = GAboundCheck(x, UB, LB)
%Boundary check, this function refelect the population if out of bound

numV = numel(x);
for k = 1:numV
    if x(k) > UB(k)
        x(k) = UB(k)-(x(k)-UB(k));
        if x(k) < LB(k)
            x(k) = UB(k);
        end
    elseif x(k) < LB(k)
        x(k) = LB(k)+(LB(k)-x(k));
        if x(k) > UB(k)
            x(k) = LB(k);
        end
    end
end


end


