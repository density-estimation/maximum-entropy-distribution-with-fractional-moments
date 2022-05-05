function [fitnessBest, popBest, fitnessCurrentBest, popCurrentBest] = ...
    de(costFun, D, sInitial, UB, LB, NP, popSize, absTol, maxGeneration, stallGeneration)
%% Input defintion
% D = variable dimension, sInitial = Initial decision vectors, ...
% UB = upper bound, LB = lower bound, NP = population size, ...
% CR = crossover probability, popSize = differential weight
% absTol = absolute tolerance
% NP = Group Size


% fitnessCurrent_DE = zeros(maxGeneration, popSize); % A matrix that stores all the fitness value
fitnessCurrent = zeros(1, popSize);  % Temp stored fitness value for each generation
sTemp = zeros(popSize,D); % Temp agent/decision vector in a certain population
fitnessCurrentBest = zeros(1, maxGeneration); % The best fitness value from each generation
popCurrentBest = zeros(maxGeneration, D); % The best decision vector/agent corresponding to the best fitness value


popCurrent = sInitial;  % Assign initial population to current population
for countNP = 1:popSize 
    % Loop from population and obtain the fitness value for 
    %each element
    fitnessCurrent(countNP) = costFun(popCurrent(countNP,:));
end



groupNumber = zeros(1,popSize);
for countPopSize = 1:popSize
    groupNumber(countPopSize) = floor((countPopSize-1)/NP)+1;
end


countG = 1; % Generation number count;
flag1 = 1; % for checking generation number 
flag2 = 1; % for checking absolute tol


while flag1&&flag2
    
    CR = 1-log(countG)/log(maxGeneration);
    
    % When generation nubmer reaches limit or absolute 
    %error is smaller than the specified tol error
    for countNP = 1:popSize  % Loop from elements in the population 
        
        
        %% Pick three agents a, b and c from population at random ...
        % They must be different from each other as well as current agent
        % countNP
        currentGroup = groupNumber(countNP);
        
        a = (currentGroup-1)*NP+randi(NP,1);
        while (a == countNP)
            a = (currentGroup-1)*NP+randi(NP,1);            
        end
        b = (currentGroup-1)*NP+randi(NP,1);
        while (b == countNP || b == a)
            b = (currentGroup-1)*NP+randi(NP,1);
        end
        c = (currentGroup-1)*NP+randi(NP,1);
        while (c == countNP || c==b || c == a)
            c = (currentGroup-1)*NP+randi(NP,1);
        end
              
        %% Pick a random index R to make sure Rth element is changed
        
        R = randi(D,1);
        F = rand(1,D)*2;
        sTemp(countNP,:) = popCurrent(a, :)+F.*(popCurrent(b,:) - popCurrent(c,:));
        
        for countD = 1:D
            if sTemp(countNP,countD) > UB(countD)
                    sTemp(countNP,countD) = UB(countD)-(sTemp(countNP,countD)-UB(countD));
                    if sTemp(countNP,countD) < LB(countD)
                        sTemp(countNP,countD) = UB(countD);
                    end
                elseif sTemp(countNP,countD) < LB(countD)
                    sTemp(countNP,countD) = LB(countD)+(LB(countD)-sTemp(countNP,countD));
                    if sTemp(countNP,countD) > UB(countD)
                         sTemp(countNP,countD) = LB(countD);
                    end
            end
                          
            if rand > CR && countD ~= R  
                sTemp(countNP,countD) = popCurrent(countNP, countD);
            end
        end
        
        
        if countNP == popSize
            for countNPT = 1:popSize
                fitnessTemp(countNPT) = costFun(sTemp(countNPT,:));
            end
        end
        
        
        if countNP == popSize
            for countNPT = 1:popSize
                if fitnessTemp(countNPT) < fitnessCurrent(countNPT)
                    fitnessCurrent(countNPT) = fitnessTemp(countNPT);
                    popCurrent(countNPT,:) = sTemp(countNPT,:);
                end  
            end
        end
        

    end
    
    
    [fitnessCurrentSort, fitnessCurrentSortIndex] = sort(fitnessCurrent);
    
    fitnessCurrent(1:NP) = fitnessCurrentSort(1:NP);
    popCurrent(1:NP,:) = popCurrent(fitnessCurrentSortIndex(1:NP),:);
    
    
%     fitnessCurrent_DE(countG,:) = fitnessCurrent;
    
    
    fitnessCurrentBest(countG) = fitnessCurrentSort(1);
    fitnessBestIndex = fitnessCurrentSortIndex(1);
    
    popCurrentBest(countG,:) = popCurrent(fitnessBestIndex,:);   
    
    fitnessBest = fitnessCurrentBest(countG);
    popBest = popCurrentBest(countG,:);
    
    flag1 = countG < maxGeneration;  
    %flag1 controls current generation number is smaller than max
    %generation number
    
    %% flag2 controls the relative error
    
    if countG <= stallGeneration
        flag2 = 1;
    else 
       flag2 = abs(fitnessBest-fitnessCurrentBest(countG-stallGeneration)) > absTol;
    end
    countG = countG +1;
   
    
end


end
    