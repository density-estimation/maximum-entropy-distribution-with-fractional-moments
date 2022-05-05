function mmodel = med_train(theta)

% med_train  Parameter estimation for Maximum Entropy probability
%            distribution with fractional moments
%
% Call:   mmodel = med_train(theta)
%
% Input
% theta  :random sample
%
% Output
% m_model:trained model with alpha, lambda, b
%
% xiaodong.zhang@u.nus.edu
% Last update December 26, 2021
% MATLAB version R2020b

%% GA algorithm parameter setting
imagineY = 100;  %the threshold which is unrealistic
stallGeneration = 200;
absTol = 1e-5;  % the different between stallGeneration number of generation
maxGeneration = 2000;  %As defined
pCrossover = 0.9;  %the probability of the occurance of Crossover
% pMutation = 0.1;
popSize = 20;   % population size for each generation
nElitism = 5; %The number of population that could be passed on to the next generation
selectionMethod = 'Tournament';  %can choose from 'Tournament', 'Roulette'

%% ---------------Parametric setting for the proposed method---------------
N = numel(theta);           % Sample size
maxNumV = 6;                % Maximum number of alpha
interUB = 100;  %Upperbound of intergration domain, Lower bound is unknonw

%% Iteration over each set of data

warning('off','all')

% Data normalization
x_DataR = theta;
x_LB = min(x_DataR);
x_UB = max(x_DataR);
x_Data = (x_DataR-x_LB)./(x_UB-x_LB);


x_Min = min(x_Data);
x_Std = std(x_Data);
x_Data(x_Data==x_Min)=[];

lambda0 = zeros(maxNumV, 1);
costBest = zeros(maxNumV, 1);
sBest_numV = zeros(maxNumV,2*maxNumV+1);
AIC = zeros(maxNumV, 1);

cost =@(x) med_cost(x, x_Data,imagineY,interUB); %Initialize cost function

for numV = 1:maxNumV
    % Parametric setting for MEGA

    vDimension = numV+1; %Number of variables
    alpha_UB = 4; %Upbound for alpha
    alpha_LB = -4; %Lowerbound for alpha
    UB = [repmat(alpha_UB,1,numV),x_Min];
    LB = [repmat(alpha_LB,1,numV),x_Min-10*x_Std];
    alpha_UB_Initial = 1; %Upbound for alpha
    alpha_LB_Initial = -1; %Lowerbound for alpha
    UB_Initial = [repmat(alpha_UB_Initial,1,numV),x_Min];
    LB_Initial = [repmat(alpha_LB_Initial,1,numV),x_Min-1*x_Std];
    
    sInitial = cat(1,repmat((UB_Initial-LB_Initial),popSize,1)).*...
        rand([popSize,vDimension])+cat(1,repmat(LB_Initial,popSize,1));
    
    [costBest(numV), alphaBest] = ...
        ga(cost, sInitial,UB,LB, popSize, maxGeneration, selectionMethod, ...
        pCrossover, nElitism, vDimension, absTol, stallGeneration, imagineY);
    
    while min(abs(alphaBest(1:end-1)-UB(1:end-1)))<1e-3 || min(abs(alphaBest(1:end-1)-LB(1:end-1)))<1e-3
        alpha_UB = alpha_UB+1; %Upbound for alpha
        alpha_LB = alpha_LB-1; %Lowerbound for alpha
        UB = [repmat(alpha_UB,1,numV),x_Min];
        LB = [repmat(alpha_LB,1,numV),x_Min-10*x_Std];
        
        sInitial = cat(1,repmat((UB-LB),popSize,1)).*...
            rand([popSize,vDimension])+cat(1,repmat(LB,popSize,1));
        [costBest(numV), alphaBest] = ...
            ga(cost, sInitial,UB,LB, popSize, maxGeneration, selectionMethod, ...
            pCrossover, nElitism, vDimension, absTol, stallGeneration, imagineY);
    end
    [~, lambdaBest,lambda0(numV)] = cost(alphaBest);
    sBest_numV(numV,1:2*numV+1) = [lambdaBest',alphaBest];
    AIC(numV) = costBest(numV)+(2*numV+1)/N; 
end

[~,costBestIndex] = min(AIC);

mmodel.interUB = interUB;
mmodel.x_LB = x_LB;
mmodel.x_UB = x_UB;
mmodel.lambda0 = lambda0(costBestIndex);
mmodel.lambda = sBest_numV(costBestIndex,1:costBestIndex);
mmodel.alpha = sBest_numV(costBestIndex,costBestIndex+1:2*costBestIndex);
mmodel.b = sBest_numV(costBestIndex,2*costBestIndex+1);

end




