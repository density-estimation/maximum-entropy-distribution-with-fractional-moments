function [y,lambda,lambda0] = med_cost(a,x_Data, imagineY, interUB) 

% med_cost :objective function
%
% Call:   [y,lambda,lambda0] = med_cost(a,x_Data, imagineY, interUB) 
%
% Input
% a:        alpha and b
% x_Data:   normalized random sample
% imagineY: an infeasible number
% interUB:  integration upper bound
%
% Output
% y:        objective function value
% lambda:   parameter
% lambda0:  normalization parameter
%
% xiaodong.zhang@u.nus.edu
% Last update August 05, 2021
% MATLAB version R2020b


lastwarn('');
lambda0 = imagineY;
y = imagineY;
numV = numel(a)-1; 

lambda = zeros(numV,1);

b = a(end);
alpha = a(1:end-1);
interLB = b;

x_Data_b = x_Data-b;

x_Data_b(x_Data_b==0)=[]; 

numData = numel(x_Data_b);

belta = zeros(numV,numV);
miu = zeros(numV,1);
%% Matrix belta and miu

for i = 1:numV
    for k = 1:numV
        belta(i,k) = alpha(k)./(alpha(i)+1).*(sum(x_Data_b.^(alpha(i)+alpha(k)))./numData);
    end
    miu(i,1) =  sum(x_Data_b.^alpha(i))./numData;
end

% try
     
if  rank(belta)<numV    %Check whether the area of pdf is zero
    return;
end

% catch
%     save('Error.mat', 'a','x_Data', 'imagineY', 'interUB');
% end

lambda= belta\miu;

sumSeries = @(x) lambda'*(bsxfun(@power, x(:), alpha))';
pdf_area = integral(@(x) exp(-sumSeries(x-b)),interLB,interUB);

[~, warnId] = lastwarn();
if  pdf_area < eps || ~(isempty(warnId))    %Check whether the area of pdf is zero
    return;
end

lambda0 = log(pdf_area); %integral from 0 to interUB

f_PDF = @(x) exp(-sumSeries(x-b))/pdf_area;
f_CDF = integral(f_PDF,interLB,interUB);

if abs(f_CDF-1) >= eps(f_CDF) || isinf(f_CDF)
    return;
end

%---------Target Minimization Function------------
y = lambda0+sum(bsxfun(@power, x_Data_b(:), alpha))./numData*lambda;

end












