clear;
clc;

mu = 10;
sigma = 1;
theta = normrnd(mu,sigma,[1,1e4]);

x = linspace(5,15,100); 
y = normpdf(x,mu,sigma);
plot(x,y)

mmodel = med_train(theta);
[y_pdf, y_cdf, y_poe] = med_predictor(mmodel, x);

% Probability density function
figure;
plot(x, normpdf(x,mu,sigma), 'DisplayName', 'Normal');
hold on 
plot(x, y_pdf, '-.', 'DisplayName', 'MED');
xlabel('x')
ylabel('Probability density')
legend();

% Cumulative distribution function 
figure;
plot(x, normcdf(x,mu,sigma), 'DisplayName', 'Normal');
hold on 
plot(x, y_cdf, '-.', 'DisplayName', 'MED');
xlabel('x')
ylabel('Cumulative distribution function')
legend();

% Probability of exceedance  
figure;
semilogy(x, normcdf(x,mu,sigma,'upper'), 'DisplayName', 'Normal');
hold on 
plot(x, y_poe, '-.', 'DisplayName', 'MED');
xlabel('x')
ylabel('Cumulative distribution function')
legend();


