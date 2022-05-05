function [pdfX, xPlot, app_name] = costf(cf_number)
% Give costf function number, return distribution name, paramter, and plot
% location

% input
% cf_number: cost function number

% output
% pdfX: the distribution parameters
% xPlot: plot points for the rancom variable


np = 100;  %number of points

switch cf_number
    case 1       % Normal distribution
        pdfX = {'Normal',-5,2};
        x_plot_lb = -16;
        x_plot_ub = 6;
        xPlot = x_plot_lb:(x_plot_ub-x_plot_lb)/(np-1):x_plot_ub;
        app_name = 'Normal';
        
    case 2       % Normal distribution
        pdfX = {'gev', 0, 1, 0};
        x_plot_lb = -5;
        x_plot_ub = 15;
        xPlot = x_plot_lb:(x_plot_ub-x_plot_lb)/(np-1):x_plot_ub;
        app_name = 'Gumbel';
        
    case 3       % Exponential distribution
        pdfX = {'exp',1};
        x_plot_lb = 0;
        x_plot_ub = 18;
        xPlot = x_plot_lb:(x_plot_ub-x_plot_lb)/(np-1):x_plot_ub;
        app_name = 'Exp';
        
    case 4       % Gamma distribution
        % k=2, b=1 (standardized)
        pdfX = {'gam',2,1};
        x_plot_lb = 0;
        x_plot_ub = 18;
        xPlot = x_plot_lb:(x_plot_ub-x_plot_lb)/(np-1):x_plot_ub;
        app_name = 'Gamma';
        
    case 5       % Inverse Gaussian
        pdfX = {'inversegaussian',1,2};
        x_plot_lb = 0;
        x_plot_ub = 16;
        xPlot = x_plot_lb:(x_plot_ub-x_plot_lb)/(np-1):x_plot_ub;
        app_name = 'InverseGaussian';
        
    case 6       % Weibull
        pdfX = {'wbl',1,2};
        x_plot_lb = 0;
        x_plot_ub = 4.5;
        xPlot = x_plot_lb:(x_plot_ub-x_plot_lb)/(np-1):x_plot_ub;
        app_name = 'Weibull';
        
    case 7       % Lognormal
        m = 1.5;
        v = 1;
        mu = log((m^2)/sqrt(v+m^2));
        sigma = sqrt(log(v/(m^2)+1));
        pdfX = {'logn',mu,sigma};
        x_plot_lb = 0;
        x_plot_ub = 30;
        xPlot = x_plot_lb:(x_plot_ub-x_plot_lb)/(np-1):x_plot_ub;
        app_name = 'Lognormal';
        
    case 8       % Chi2
        v = 15;
        pdfX = {'chi2',v};
        x_plot_lb = 0;
        x_plot_ub = 70;
        xPlot = x_plot_lb:(x_plot_ub-x_plot_lb)/(np-1):x_plot_ub;
        app_name = 'Chi2';
        
    case 9       % F
        pdfX = {'f',5,15};
        x_plot_lb = 0;
        x_plot_ub = 18;
        xPlot = x_plot_lb:(x_plot_ub-x_plot_lb)/(np-1):x_plot_ub;
        app_name = 'F';
        
        
    case 10       % T
        pdfX = {'T',5};
        x_plot_lb = -25;
        x_plot_ub = 25;
        xPlot = x_plot_lb:(x_plot_ub-x_plot_lb)/(np-1):x_plot_ub;
        app_name = 'T';
        
    otherwise
        disp('other value')
end
