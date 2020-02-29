%% vigsiv-CSS-L-STOC-ECF: Figure 2: Convergence
% This code runs the smoothing example in Fig. 2 presented in the paper: 
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
plot_fontSize = 9;
plot_linewidth = 2;

% Confidence interval params: 

alpha = 0.2;

%% Construct data vector:
nt = 1E5;
datatrue = [unifrnd(0,5,nt,1); gamrnd(2,5,[nt,1])]';
n1 = 10;
data1 = [unifrnd(0,5,n1,1); gamrnd(2,5,[n1,1])]';
n2 = 100;
data2 = [unifrnd(0,5,n2,1); gamrnd(2,5,[n2,1])]';
n3 = 1000;
data3 = [unifrnd(0,5,n3,1); gamrnd(2,5,[n3,1])]';

%% Compute True CDF, i.e. ECF with large number of samples: 
[sigmat,~,~,~] = kde(datatrue,nt,min(datatrue),max(datatrue));
cf_func = @(t) diracMixture(t,datatrue,sigmat);
clear options
options.isPlot = false;
result_true = cf2DistGP(cf_func,[],[],options);

%% Compute True CDF, i.e. ECF with n = 10: 
[sigma1,~,~,~] = kde(data1,n1,min(data1),max(data1));
cf_func = @(t) diracMixture(t,data1,sigma1);
clear options
options.isPlot = false;
result1 = cf2DistGP(cf_func,[],[],options);
conf_int1 = sqrt(log(2/alpha)/(2*n1));

%% Compute True CDF, i.e. ECF with n = 100: 
[sigma2,~,~,~] = kde(data2,n2,min(data2),max(data2));
cf_func = @(t) diracMixture(t,data2,sigma2);
clear options
options.isPlot = false;
result2 = cf2DistGP(cf_func,[],[],options);
conf_int2 = sqrt(log(2/alpha)/(2*n2));


%% Compute True CDF, i.e. ECF with n = 1000: 
[sigma3,~,~,~] = kde(data3,n3,min(data3),max(data3));
cf_func = @(t) diracMixture(t,data3,sigma3);
clear options
options.isPlot = false;
result3 = cf2DistGP(cf_func,[],[],options);
conf_int3 = sqrt(log(2/alpha)/(2*n3));

%% Compute Moments: 

H = [1000:1E3:1E6];

for j = 1:length(H)
        tt = tic;
        datam = [];
        datam = [unifrnd(0,5,H(j),1); gamrnd(2,5,[H(j),1])]';
        m(j) = diracMixtureCostmean(datam);
        [sigma,~,~,~] = kde(datam,H(j),min(datam),max(datam));
        m2(j) = diracMixtureCostcov(datam,sigma);
        ttf(j) = toc(tt);
end

stdm = std(m,1,2);
stdm2 = std(m2,1,2);



%% Plotting: 

Fig2 = figure('Units', 'points', ...
       'Position', [0, 0, width, height]);
ax = axes;
ax.Units = 'points';

subplot(2,3,1)
histogram(datatrue,'Normalization','cdf');
hold on
plot(result_true.x,result_true.cdf,'r','LineWidth',1.5)
plot(result1.x,result1.cdf,'LineWidth',1.5)
plot(result1.x,result1.cdf+conf_int1,'-b','LineWidth',1.5)
plot(result1.x,result1.cdf-conf_int1,'-b','LineWidth',1.5)
% ax = gca;
% ax.Position = [0.1 0.1 0.4 0.8];
axis([-13 20 0 1])
xlabel('$x$')
ylabel('$\Phi_{\textbf{w}}(x)$')

subplot(2,3,2)
histogram(datatrue,'Normalization','cdf');
hold on
plot(result_true.x,result_true.cdf,'r','LineWidth',1.5)
plot(result2.x,result2.cdf,'LineWidth',1.5)
plot(result2.x,result2.cdf+conf_int2,'-b','LineWidth',1.5)
plot(result2.x,result2.cdf-conf_int2,'-b','LineWidth',1.5)
% ax = gca;
% ax.Position = [0.1 0.1 0.4 0.8];
axis([-13 20 0 1])
xlabel('$x$')
set(gca,'ytick',[])

subplot(2,3,3)
histogram(datatrue,'Normalization','cdf');
hold on
plot(result_true.x,result_true.cdf,'r','LineWidth',1.5)
plot(result3.x,result3.cdf,'LineWidth',1.5)
plot(result3.x,result3.cdf+conf_int3,'-b','LineWidth',1.5)
plot(result3.x,result3.cdf-conf_int3,'-b','LineWidth',1.5)
% ax = gca;
% ax.Position = [0.1 0.1 0.4 0.8];
axis([-5 25 0 1])
xlabel('$x$')
set(gca,'ytick',[])

subplot(2,3,4)
plot(H,m(:,1)) %,stdm)
hold on
yline(m(end,1),'-r','Linewidth',1.5);
ylabel('E[\textbf{w}]')
xlabel('\# samples')
set(gca,'xscale','log')

subplot(2,3,5)
plot(H,m2(:,1)) %stdm2)
hold on
yline(m2(end,1),'-r','Linewidth',1.5);
ylabel('E[\textbf{w}$^2$]')
xlabel('\# samples')
set(gca,'xscale','log')



function cf = diracMixture(t,data,sigma)

	t = reshape(t,length(t),1);
    cf_int = sum(1/size(data,2)*exp(1i *  t * data),2).*exp(-(sigma*t).^2/2);
    cf = cf_int; 

end

function m = diracMixtureCostmean(data)

	t = 0;

    m = (1i)^(-1)*1/(size(data,2))*sum(1i*data,2);


end

function m2 = diracMixtureCostcov(data,sigma)

	t = 0;

    m2 = (1i)^(-2)*1/(size(data,2))*...
        sum(-sigma.^2+(1i*data).^2,2);

end