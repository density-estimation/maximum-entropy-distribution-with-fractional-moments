function [y_pdf, y_cdf, y_poe] = med_predictor(mmodel, x)

% med_predictor  Probability density estimation from maximum entropy
%                distribution with fractional moments model
%
% Call:   [y_pdf, y_cdf, y_poe] = med_predictor(mmodel, x)
%
% Input
% mmodel  : a structure containing model parameters
%
% Output
% y_pdf : probability density
% y_cdf : cumulative density
% y_poe : probability of exceedance
%
% xiaodong.zhang@u.nus.edu
% Last update August 05, 2021
% MATLAB version R2020b

xLB = mmodel.xLB;
xUB = mmodel.xUB;
lambda0 = mmodel.lambda0;
lambda = mmodel.lambda;
alpha = mmodel.alpha;
b = mmodel.b;
interUB = mmodel.interUB;

nP = numel(x);

y_pdf = zeros(1,nP);
y_cdf = zeros(1,nP);
y_poe = zeros(1,nP);
    
sumSeries = @(x) lambda*(bsxfun(@power, x(:), alpha))';

f_PDF = @(x) exp(-(lambda0+sumSeries(x-b)));
f_CDF = @(x) integral(f_PDF,b,x);
f_POE = @(x) integral(f_PDF,x,interUB);

xN = (x-xLB)./(xUB-xLB);
for k = 1:nP
    y_pdf(k) = f_PDF(xN(k))./(xUB-xLB);
    y_cdf(k) = f_CDF(xN(k));
    y_poe(k) = f_POE(xN(k));
end

y_pdf = real(y_pdf);
y_cdf = real(y_cdf);
y_poe = real(y_poe);

end