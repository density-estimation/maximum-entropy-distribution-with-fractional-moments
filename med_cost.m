function [y,lambda,lambda0] = med_cost(a,x_Data, infY, interUB)

% med_cost :objective function
%
% Call:   [y,lambda,lambda0] = med_cost(a,x_Data, infY, interUB)
%
% Input
% a:        alpha and b
% x_Data:   normalized random sample
% infY:     an infeasible number
% interUB:  integration upper bound
%
% Output
% y:        objective function value
% lambda:   parameter
% lambda0:  normalization parameter
%
% xiaodong.zhang@u.nus.edu
% Last update Jan 1, 2022
% MATLAB version R2020b


lastwarn('');
lambda0 = infY;
y = infY;
numV = numel(a)-1;

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

% if  rank(belta)<numV    %Check whether the area of pdf is zero
%     return;
% end

lambda= belta\miu;


% %% Newton's method
% xLamda0 = lambda;
% 
% % Matrix H   
% 
% tolAbs = 1e-3;      %tolerance for relative error
% tolRel = 1e-3;
% errRel = 1;
% errAbs = 1;
% flag1 = errAbs>tolAbs && errRel>tolRel; %flag1 is for relative error
% 
% miuCount = 0;         
% 
% while flag1 && miuCount<15 
%     miuCount = miuCount+1;
%     sumSeries = @(x) xLamda0'*(bsxfun(@power, x(:), alpha))';
%     pdf_area = integral(@(x) exp(-sumSeries(x-b)),interLB,interUB);
%     
%     % Matrix miuIJ and miuIpJ
%     miuIpJ = zeros(numV, numV);
%     miuLamda = zeros(numV,1);
%     for i = 1:numV
%         for j = i:numV
%             fun = @(x) x.^(alpha(i)+alpha(j)).*exp(-sumSeries(x-b));
%             miuIpJ(i,j) = integral(fun,interLB,interUB)./pdf_area;  %first row
%             miuIpJ(j,i) = miuIpJ(i,j);
%         end      
%         miuLamda(i,1) = (integral(@(x) x.^(alpha(i)).*exp(-sumSeries(x-b)),...
%             interLB,interUB))./pdf_area; 
%     end
%                                        
%     Gradiant = miu-miuLamda;
%     %% Hessian maatrix should be positive definite
%     miuItJ = zeros(numV, numV);
%     for i = 1:numV
%         for j = i:numV
%             miuItJ(i,j) = miuLamda(i).* miuLamda(j);
%             miuItJ(j,i) = miuItJ(i,j);
%         end
%     end 
%    
%     %% Hessian maatrix should be positive definite
%  
%     H = miuIpJ-miuItJ;
%     [~, cholp] = chol(H);
%     if cholp~=0
%         break;
%     end
% 
%     xLamda1 = xLamda0-H\Gradiant;
%     errAbs = norm(xLamda1-xLamda0,inf);
%     errRel = norm((xLamda1-xLamda0)./xLamda0,inf);
%     xLamda0 = xLamda1;
%     
%     flag1 = errAbs>tolAbs && errRel>tolRel;
% 
% end
% 
% % Cost function value
% 
% if  cholp~=0 || flag1==1  %Hessian matrix is not positive definite
%     lambda =lambda;
% else
%     lambda = xLamda0;
% end


sumSeries = @(x) lambda'*(bsxfun(@power, x(:), alpha))';
pdf_area = integral(@(x) exp(-sumSeries(x-b)),interLB,interUB);

if isinf(pdf_area) || (pdf_area < eps)
    return;
end

lambda0 = log(pdf_area); %integral from 0 to interUB

f_PDF = @(x) exp(-sumSeries(x-b))/pdf_area;
f_CDF = integral(f_PDF,interLB,interUB);

[~, warnId] = lastwarn();

%---------Target Minimization Function------------
y_t = lambda0+sum(bsxfun(@power, x_Data_b(:), alpha))./numData*lambda;

if  abs(f_CDF-1) > 5*eps || ~(isempty(warnId))
    y = 0.9*y_t;
else
    y = y_t;
end


end


