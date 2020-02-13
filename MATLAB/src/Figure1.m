%% vigsiv-CSS-L-STOC-ECF: Figure 1: Smoothing
% This code runs the smoothing example in Fig. 1 presented in the paper: 
% "Stochastic Optimal Using Empirical Characteristic Functions"
%
% REQUIRED DEPENDENCIES: - CharFunTool 
%                          (https://github.com/witkovsky/CharFunTool/)
%                        - MATLAB Statistics and Machine Learning
%                          Toolbox


%% Housekeeping 
clc, clear, close all


% Figure params: 

width = 252; 
height = 100;
plot_markersize = 15;
plot_fontSize = 8;
plot_linewidth = 2;

%% Construct data vector:
rng(101)
n = 300;
data = [unifrnd(-4,5,n,1); exprnd(2,n,1)]';

%% Compute CDF when there is no smoothing, \sigma = 0
sigma = 0;
cf_func = @(t) diracMixturens(t,data);
clear options
options.isPlot = false;
resultns = cf2DistGP(cf_func,[],[],options);


%% Compute PDF/CDF when we algorithmically chosen, \sigma = algorithm
[sigma,~,~,~] = kde(data,n,min(data),max(data));
cf_func = @(t) diracMixture(t,data,sigma);
clear options
options.isPlot = false;
result = cf2DistGP(cf_func,[],[],options);
fig2d = figure('Units', 'points', ...
       'Position', [0, 0, width, height]);
ax = axes;
ax.Units = 'points';

subplot(2,2,[1 3])
histogram(data,'Normalization','cdf');
hold on
plot(resultns.x,resultns.cdf,'r','LineWidth',1.5)
% ax = gca;
% ax.Position = [0.1 0.1 0.4 0.8];
axis([-6 13 0 1])
xlabel('$x$')
ylabel('$\Phi_{\textbf{w}}(x)$')

subplot(2,2,2)
histogram(data,'Normalization','cdf');
hold on
plot(resultns.x,resultns.cdf,'r','LineWidth',1.5)
set(gca,'xtick',[])
set(gca,'ytick',[])
box on
axis([7.5 11.5 0.986 1])

subplot(2,2,4)
set(findall(gcf,'-property','FontSize'),'FontSize',plot_fontSize)
histogram(data,'Normalization','cdf');
hold on
plot(result.x,result.cdf,'r','LineWidth',1.5)
set(gca,'xtick',[])
set(gca,'ytick',[])
box on
axis([7.5 11.5 0.986 1])


function cf = diracMixture(t,data,sigma)

	t = reshape(t,length(t),1);
    cf_int = sum(1/size(data,2)*exp(1i *  t * data),2).*exp(-(sigma*t).^2/2);
    cf = cf_int; 

end

function cf = diracMixturens(t,data)

	t = reshape(t,length(t),1);
    cf_int = sum(1/size(data,2)*exp(1i *  t * data),2);
    cf = cf_int; 

end