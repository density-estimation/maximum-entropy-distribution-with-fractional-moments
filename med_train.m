function mmodel = med_train(theta,b)

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
% Last update Jan 1, 2022
% MATLAB version R2020b

%% ---------------Parametric setting for the proposed method---------------
warning('off','all')
optM = 'GeneticAlgorithm';
% optM = 'DifferentialEvolution';

N = numel(theta)-1;                                                        % Sample size
maxM = 8;                                                                  % Maximum number of alpha
interUB = 10;                                                              % Upperbound of intergration domain, Lower bound is unknonw
% Data normalization
xLB = min(theta);
xUB = max(theta);
xNorm = (theta-xLB)./(xUB-xLB);                                            % Normalized data
xStd = std(xNorm);
xNorm(xNorm==0)=[];                                                        % Required, not redundant


infY = 100;
% Optimization parameter
absTol = 1./N*1e-2;                                                        % difference between stallGeneration number of generations
stallGeneration = 100;                                                      % number of stall generation
maxGeneration = 2000;                                                      % Maximum number of generations
popSize = 20;

if strcmp('u',b)
    
    cost =@(x) med_cost(x,xNorm,infY,interUB); %Initialize cost function
    
    switch optM
        case 'GeneticAlgorithm'
            %%-------------------GA algorithm parameter setting------------------------
            pCrossover = 0.9;                                                          % Crossover probability
            nElitism = 5;                                                              % Number of population that could be passed on to the next generation
            selectionMethod = 'Tournament';                                            % 'Tournament', 'Roulette' selection method
            opt = @(vDimension, sInitial, UB, LB) ga(cost, vDimension, sInitial,...
                UB,LB,popSize, maxGeneration, selectionMethod, pCrossover, nElitism,...
                absTol, stallGeneration);
        case 'DifferentialEvolution'
            NP = 5;
            opt = @(D, sInitial, UB, LB) de(cost, D, sInitial, UB, LB, NP, popSize,...
                absTol, maxGeneration, stallGeneration);
    end
    
    % Iteration over each set of data
    costBest = zeros(maxM, 1);
    AIC = zeros(maxM, 1);
    sB = cell(maxM, 1);
    
    for m = 1:maxM
               
        d = m+1;                                                               % Number of variables
        alphaUB = 4;                                                           % Upperbound for alpha
        alphaLB = -4;                                                          % Lowerbound for alph
        t = 1;
        bUB = 0;
        bLB = -t*xStd;
        UB = [repmat(alphaUB,1,m),bUB];
        LB = [repmat(alphaLB,1,m),bLB];
        
        sInitial = cat(1,repmat((UB-LB),popSize,1)).*lhsdesign(popSize,d)+...
            cat(1,repmat(LB,popSize,1));
        [costBest(m), sBest] = opt(d, sInitial, UB, LB);
        
        abstol = 1e-3;
        flag1 = min(abs(sBest(1:m)-alphaUB)) < abstol;
        flag2 = min(abs(sBest(1:m)-alphaLB)) < abstol;
        flag3 = min(abs(sBest(m+1)-bLB)) < abs(bLB*abstol);
        
        while  flag1 || flag2 || flag3
            
            if flag1 || flag2
                alphaUB = alphaUB+1;                                           % Upperbound for alpha
                alphaLB = alphaLB-1;                                           % Lowerbound for alpha
            else
                t = t+1;
                bLB = -t*xStd;
            end
            UB = [repmat(alphaUB,1,m),bUB];
            LB = [repmat(alphaLB,1,m),bLB];
            
            sInitial = cat(1,repmat((UB-LB),popSize,1)).*...
                lhsdesign(popSize,d)+cat(1,repmat(LB,popSize,1));
            [costBest(m), sBest] = opt(d, sInitial, UB, LB);
            
            flag1 = min(abs(sBest(1:m)-alphaUB)) < abstol;
            flag2 = min(abs(sBest(1:m)-alphaLB)) < abstol;
            flag3 = min(abs(sBest(m+1)-bLB)) < abs(bLB*abstol);
            
            if alphaUB >= 6 || alphaLB <=-6 || bLB<-0.6
                break
            end
        end
        sB{m,1} = sBest;
        AIC(m) = costBest(m)+(2*m+1)/N;
        if m>1 && AIC(m)>AIC(m-1)
            M = m-1;
            break;
        end
        
    end
    
    s = sB{M,1};
    [~, lambda,lambda0] = cost(s);
    alpha = s(1:M);
    b = s(M+1);
    
else
    
    b = (b-xLB)./(xUB-xLB);
    cost =@(x) med_cost([x,b],xNorm,infY,interUB); %Initialize cost function
    
    switch optM
        case 'GeneticAlgorithm'
            %%-------------------GA algorithm parameter setting------------------------
            pCrossover = 0.9;                                                          % Crossover probability
            nElitism = 5;                                                              % Number of population that could be passed on to the next generation
            selectionMethod = 'Tournament';                                            % 'Tournament', 'Roulette' selection method
            opt = @(vDimension, sInitial, UB, LB) ga(cost, vDimension, sInitial,...
                UB,LB,popSize, maxGeneration, selectionMethod, pCrossover, nElitism,...
                absTol, stallGeneration);
        case 'DifferentialEvolution'
            NP = 5;
            opt = @(D, sInitial, UB, LB) DE(cost, D, sInitial, UB, LB, NP, popSize,...
                absTol, maxGeneration, stallGeneration);
    end
    
    % Iteration over each set of data
    costBest = zeros(maxM, 1);
    AIC = zeros(maxM, 1);
    sB = cell(maxM, 1);
    
    for m = 1:maxM
        d = m;                                                               % Number of variables
        alphaUB = 4;                                                           % Upperbound for alpha
        alphaLB = -4;                                                          % Lowerbound for alpha
        UB = repmat(alphaUB,1,m);
        LB = repmat(alphaLB,1,m);
        
        sInitial = cat(1,repmat((UB-LB),popSize,1)).*lhsdesign(popSize,d)+...
            cat(1,repmat(LB,popSize,1));
        [costBest(m), sBest] = opt(d, sInitial, UB, LB);
        
        abstol = 1e-3;
        flag1 = min(abs(sBest(1:m)-alphaUB)) < abstol;
        flag2 = min(abs(sBest(1:m)-alphaLB)) < abstol;
        
        while  flag1 || flag2
            
            alphaUB = alphaUB+1;                                           % Upperbound for alpha
            alphaLB = alphaLB-1;                                           % Lowerbound for alpha
            UB = repmat(alphaUB,1,m);
            LB = repmat(alphaLB,1,m);
            
            sInitial = cat(1,repmat((UB-LB),popSize,1)).*...
                lhsdesign(popSize,d)+cat(1,repmat(LB,popSize,1));
            [costBest(m), sBest] = opt(d, sInitial, UB, LB);
            
            flag1 = min(abs(sBest(1:m)-alphaUB)) < abstol;
            flag2 = min(abs(sBest(1:m)-alphaLB)) < abstol;
            
            if alphaUB >= 6 || alphaLB <=-6
                break;
            end
        end
        sB{m,1} = sBest;
        AIC(m) = costBest(m)+(2*m+1)/N;
        if m>1 && AIC(m)>AIC(m-1)
            M = m-1;
            break;
        end
        
    end
    
    s = sB{M,1};
    [~, lambda,lambda0] = cost(s);
    alpha = s;
    
end


mmodel.interUB = interUB;
mmodel.xLB = xLB;
mmodel.xUB = xUB;
mmodel.lambda0 = lambda0;
mmodel.lambda = lambda';
mmodel.alpha = alpha;
mmodel.b = b;

end
















