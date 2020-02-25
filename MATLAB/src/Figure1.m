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

width = 505.89; 
height = 50;
plot_markersize = 15;
plot_fontSize = 8;
plot_linewidth = 2;

%% Construct data vector:
n = 1000;
data = [normrnd(0,5,n,1); wblrnd(4,2,n,1)]';
errordes = 1E-3;

%% Compute PDF/CDF when we algorithmically chosen, \sigma = algorithm
[sigma,~,~,~] = kde(data,n,min(data),max(data));
cf_func = @(t) diracMixture(t,data,sigma);
clear options
options.isPlot = false;
options.xN = 3000;
result = cf2DistGP(cf_func,[],[],options);
x = fliplr(result.x)';
cdf = fliplr(result.cdf)';

[pu_m,pu_c,res] =...
    piecewiseUnder(x,cdf,errordes,20);

pu_m = [pu_m 0];
pu_c = [pu_c cdf(end)];
xind = find(x==res(1));
y = min(pu_m.*x(xind:end)+pu_c,[],2);





fig2d = figure('Units', 'points', ...
       'Position', [0, 0, width, height]);
ax = axes;
ax.Units = 'points';

subplot(1,5,1)
histogram(data,'Normalization','cdf');
hold on
plot(result.x,result.cdf,'r','LineWidth',1.5)
% ax = gca;
% ax.Position = [0.1 0.1 0.4 0.8];
% axis([-6 13 0 1])
xlabel('$x$')
axis([-15 18 0 1])
ylabel('$\Phi_{\textbf{w}}(x)$')

subplot(1,5,2)
histogram(data,'Normalization','cdf');
hold on
plot(result.x,result.cdf,'r','LineWidth',1.5)
% ax = gca;
% ax.Position = [0.1 0.1 0.4 0.8];
% axis([-6 13 0 1])
axis([-15 18 0 1])
xlabel('$x$')
set(gca,'ytick',[])

subplot(1,5,3)
histogram(data,'Normalization','cdf');
hold on
plot(result.x,result.cdf,'r','LineWidth',1.5)
% ax = gca;
% ax.Position = [0.1 0.1 0.4 0.8];
% axis([-6 13 0 1])
axis([-15 18 0 1])
xlabel('$x$')
set(gca,'ytick',[])

subplot(1,5,4)
histogram(data,'Normalization','cdf');
hold on
plot(result.x,result.cdf,'r','LineWidth',1.5)
plot(x(xind:end),y,'g','LineWidth',1.5)
% ax = gca;
% ax.Position = [0.1 0.1 0.4 0.8];
% axis([-6 13 0 1])
axis([-15 18 0 1])
xlabel('$x$')
set(gca,'ytick',[])

subplot(1,5,5)
yline(1E-3,'-r','Linewidth',2)
hold on
plot(x(xind:end),cdf(xind:end)-y,'g','LineWidth',1.5)
axis([5 25 0 1E-3])
xlabel('$x$')
ylabel('Error')

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