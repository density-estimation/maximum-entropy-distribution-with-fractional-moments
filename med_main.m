%% Sample data generation
% Cost function number with name:
% 1, Normal 2, Gumbel 3, Exponential 4, Gamma 5, Inverse Gaussian 6, Weibull 7, Lognormal 8, Chi2
clear;
clc;
cf_number = 1;                                       % Cost function number
numRuns = 1;                                       % Nubmer of set of data
N = 2000;                                            % Sample size
[pdfX, ~, app_name] = costf(cf_number);
theta = random(pdfX{1},pdfX{2:end},[numRuns,N]);     % Only applies for theoretical distribution
folder_name = sprintf('%s/Results/%s',pwd,app_name);
if ~exist(folder_name, 'dir')
    mkdir(folder_name)
end

b = 'u';
theta_dir=sprintf('%s/%s_%d.mat',folder_name,app_name,N); %theta directory
save(theta_dir,'theta');

%% Parameter estimation
% clear;
% clc;
% cf_number = 7;
% numRuns = 6;
% N = 5000;
algorithm_name = 'med_train';
% Data set
[~, ~, app_name] = costf(cf_number);
theta_dir=sprintf('%s/Results/%s/%s_%d.mat',pwd, app_name, app_name,N);
theta = importdata(theta_dir);


mmodels = cell(numRuns, 1);
% https://www.mathworks.com/help/parallel-computing/parfor.html#d117e54445
if numRuns>1
    arg = Inf;  % Maximum number of workers
else
    arg = 0;    % Serial instead of parallel
end
tic
parfor (run_number = 1:numRuns,arg)
    fh = str2func(algorithm_name);
    mmodels{run_number}= fh(theta(run_number,:),b);
end
toc
model_dir=sprintf('%s/Results/%s/%s_%d_models.mat',pwd, app_name, app_name,N);
save(model_dir,'mmodels');

%% Results
% clear;
% clc;
% 
% cf_number = 7;
% numRuns = 6;
% N = 5000;
[~, xPlot, app_name] = costf(cf_number);

model_dir=sprintf('%s/Results/%s/%s_%d_models.mat',pwd, app_name, app_name,N);
mmodels = importdata(model_dir);

np = numel(xPlot);
y_pdf_p_a = zeros(numRuns, np);
y_cdf_p_a = zeros(numRuns, np);
y_poe_p_a = zeros(numRuns, np);

for k = 1:numRuns
    [y_pdf_p_a(k,:), y_cdf_p_a(k,:), y_poe_p_a(k,:)] = med_predictor(mmodels{k,1}, xPlot);
end

if numRuns>1
    y_pdf_p = mean(y_pdf_p_a);
    y_cdf_p = mean(y_cdf_p_a);
    y_poe_cov = std(y_poe_p_a)./mean(y_poe_p_a);
    y_poe_p = mean(y_poe_p_a);
else
    y_pdf_p = y_pdf_p_a;
    y_cdf_p = y_cdf_p_a;
    y_poe_p = y_poe_p_a;
    y_poe_cov = 0;
end

result_dir=sprintf('%s/Results/%s/%s_%d_results.mat',pwd, app_name, app_name,N);
save(result_dir,'xPlot','y_pdf_p','y_cdf_p','y_poe_p','y_poe_cov');

%% Plot
% clear;
% clc;
% clf;
% cf_number = 7;
% N = 5000;
[pdfX, xPlot, app_name] = costf(cf_number);
result_dir=sprintf('%s/Results/%s/%s_%d_results.mat',pwd, app_name, app_name, N);

figure;
y_pdf_t = pdf(pdfX{1},xPlot,pdfX{2:end});
semilogy(xPlot,y_pdf_t);
hold on
y_pdf_p = load(result_dir,'y_pdf_p').y_pdf_p;
semilogy(xPlot,y_pdf_p);

figure;
y_cdf_t = cdf(pdfX{1},xPlot,pdfX{2:end});
semilogy(xPlot,y_cdf_t);
hold on
y_cdf_p = load(result_dir,'y_cdf_p').y_cdf_p;
semilogy(xPlot,y_cdf_p);

figure;
y_poe_t = cdf(pdfX{1},xPlot,pdfX{2:end},'upper');
semilogy(xPlot,y_poe_t);
hold on
y_poe_p = load(result_dir,'y_poe_p').y_poe_p;
semilogy(xPlot,y_poe_p);


