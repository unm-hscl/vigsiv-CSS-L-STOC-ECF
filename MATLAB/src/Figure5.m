%% vigsiv-CSS-L-STOC-ECF: Figure 4: Hypersonic Vehicle
% This code runs a simple double-integrator, LTI system with UNKNOWN 
% additive disturbance to evaluate the probability that the state will lie
% in a specified envelope using approximate cdfs from the 
% characteristic function (using CharFunTool).
%
% REQUIRED DEPENDENCIES: - CharFunTool 
%                          (https://github.com/witkovsky/CharFunTool/)
%                        - SReachTools
%                          (https://unm-hscl.github.io/SReachTools/)
%                        - MATLAB Statistics and Machine Learning
%                          Toolbox

%% Housekeeping
clc, clear, close all

% Figure params: 

width = 252; 
height = 150;
plot_markersize = 15;
plot_fontSize = 8;
plot_linewidth = 2;


prob.Ts = 0.25; % Step size
prob.T = 10; % Time horizon


% Hypersonics_linearization

% A = prob.Ts*A+eye(5);
% B = prob.Ts*B;

A = [0  0        -7702.0808 7702.0808 0;
     0 -0.001346  21.284    -32.2     0;
     0 -9.919E-7  -0.0696    0        1;
     0  0         0          0        1;
     0 -2.08E-6   2.945      0        0];
 
B = [0          0;
     24.43 -37.19    ;
     -8.385E-5 -0.01121; 
     0         0;
     .123 -1.49];
 
     
sysC = ss(A,B,zeros(1,5),0); % setup a state-space model
sysD = c2d(sysC,prob.Ts); % discretize the dynamics
Adis = sysD.A;
Bdis = sysD.B; 
     
% Distrubance matrix (12x3 matrix): 
% 
% Note that there is only a disturbance on the states for the position
% of the quadcopter. 

    G = [1 0;0 1; 0 0; 0 0; 0 0];
     
     
 % Defining the disturbance: 
    
    
    dist = RandomVector('UserDefined',...
                        @(N) [zeros(2,N)]);
     
                         
                   
    sys_lti = LtiSystem('StateMatrix', Adis,...,
                'InputMatrix',Bdis,...
                'InputSpace', Polyhedron('lb',[-ones(2,1)],...
                'ub', [ones(2,1)]),...
                'DisturbanceMatrix', G,...
                'Disturbance', dist);
            
            
    
% Generate concat matrices 

    [prob.Ad, prob.Bd, prob.Gd] = getConcatMats(sys_lti, prob.T); 
    
    
% Generate data vector: 
    
n = 1000;
for i = 1:size(prob.Gd(:,1:2:end),2)
    data(2*i-1:2*i,:) = [2*wblrnd(5,4,[n,1])'; gamrnd(5,1,n,1)'];
    for j = 2*i-1:2*i
    [sigma(j),~,~,~] = kde(data(j,:),n,min(data(j,:)),max(data(j,:)));
    end
end

prob.N = 50;
for i = 1:size(prob.Gd(:,1:2:end),2)
    prob.realizations(2*i-1:2*i,:) = [2*wblrnd(5,4,[prob.N,1])'; gamrnd(5,1,prob.N,1)'];

end


% Compute upper and lower bounds on the mean: 

for i = 1:size(prob.Gd,1)
    
   
    cf_func_c = @(t) diracMixtureICC(t,prob.Gd(i,:)*data,prob.Gd(i,:)*diag(sigma)*prob.Gd(i,:)');
    clear options
    options.isPlot = false;
    options.xN = n;
    options.xMin = min(prob.Gd(i,:)*data); 
    options.xMax = max(prob.Gd(i,:)*data);
    result_conf{i} = cf2DistGP(cf_func_c,[],[],options);

    xc{i} = fliplr(result_conf{1,i}.x)';
    cdfc{i} = fliplr(result_conf{1,i}.cdf)';
    
    % Compute confidence: 
    
    confidence = 0.90;
    epsil(i) = sqrt(1/(2*n)*log(1/(1-confidence)));
    
    cdfl{i} = cdfc{i} - epsil(i); 
    cdfu{i} = cdfc{i} + epsil(i);
    
    % Mean confidence intervals: 
    
    GdWu(i) = max(xc{i}) - trapz(xc{i},cdfl{i}); 
    GdWl(i) = max(xc{i}) - trapz(xc{i},cdfu{i}); 
    
%     figure
%     title('CDF')
%     hold on
%     empcdf = histogram(prob.Gd(i,:)*data,'Normalization','cdf');
%     hold on
%     plot(xc{i},cdfc{i},'-b','LineWidth',2)
%     plot(xc{i},cdfl{i},'-r','Linewidth',2)
%     plot(xc{i},cdfu{i},'-r','Linewidth',2)
    
end
    
    
% Compute the cdf using the Gil-Pelaez inversion formula: 

    % Specify safe set: 

    p = [-1 zeros(1,4);
          1 zeros(1,4);
          0 -1  0 0 0;
          0  1  0 0 0;];
    q = [-85000 85200 -7650 7750]';
    
%     p = [-1  0 0   0 0;
%           1  -1 0   0 0;
%           0  1 -1 0 0;
%           0  0  1 0 0;];
%     q = [84800 85200 deg2rad(10) deg2rad(10)]';

    pbig = kron(eye(prob.T),p);
    qbig = repmat(q,prob.T,1);

    prob.pbig = pbig; 
    prob.qbig = qbig;

    
    n_lin_const = size(pbig,1); 


disp('Generating PWAUA for ICC')    
tic    
for k = 1:n_lin_const
                
                transformed_rv = prob.pbig(k,:)*prob.Gd*data;
                
                cf_func = @(t) diracMixtureICC(t,transformed_rv,(prob.pbig(k,:)*prob.Gd)*diag(sigma)*(prob.pbig(k,:)*prob.Gd)');
                clear options
                options.isPlot = false;
                options.xN = 1000; 
                options.xMin = min(transformed_rv); 
                options.xMax = max(transformed_rv);
%                 options.N = 1000;
                result{k} = cf2DistGP(cf_func,[],[],options);
                
                x{k} = fliplr(result{1,k}.x)';
                pdf{k} = fliplr(result{1,k}.pdf)';
                cdf{k} = fliplr(result{1,k}.cdf)';
                
%                 x{k} = result{1,k}.x';
%                 pdf{k} = result{1,k}.pdf';
%                 cdf{k} = result{1,k}.cdf';

                [prob.pu_m{k},prob.pu_c{k},res(k,:)] =...
                    piecewiseUnder(x{k},cdf{k},1E-3,20);
                
                prob.pu_m{k} = [prob.pu_m{k} 0];
                prob.pu_c{k} = [prob.pu_c{k} 1];
                
                xind = find(x{k}==res(k,1));
                
                y{k} = min(prob.pu_m{k}.*x{k}(xind:end)+prob.pu_c{k},[],2);
                    
%                 figure
%                 title('CDF')
%                 hold on
%                 empcdf = histogram(transformed_rv,'Normalization','cdf');
%                 hold on
%                 plot(x{k},cdf{k},'-b','LineWidth',2)
%                 plot(x{k}(xind:end),y{k},'-r','Linewidth',2)
%                 figure(2); 
%                 plot(x{k}(xind:end),cdf{k}(xind:end)-y{k})
%                 hold on
                
end
toc

%% Problem params: 

prob.Delta = 1-0.8;
prob.x0 = [85000; 7702.0808;deg2rad(1.5153);deg2rad(1.5153);0;];
Xd1 = linspace(85000,85000,prob.T)';
Xd2 = linspace(7702.0808,7702.0808,prob.T)';
Xd3 = linspace(deg2rad(1.5153),deg2rad(1.5153),prob.T)';
Xd4 = linspace(deg2rad(1.5153),deg2rad(1.5153),prob.T)';
Xd5 = linspace(0,0,prob.T)';
prob.Xd = reshape([Xd1'; Xd2';Xd3';Xd4';Xd5'],[],1);
prob.ulimu = repmat([1.2; 0.2618],prob.T,1);
prob.uliml = repmat([0.2; -0.2618],prob.T,1);
prob.xlb = res(:,1);

% prob.xterm = prob.x0;

%% Generate Moments for the Quadratic Cost: 
prob.R = 1E-2*eye(size(prob.Bd,2));
prob.Q = 10*eye(size(prob.Ad,1));
% prob.Q(end-4:end,end-4:end) = 1000*prob.Q(end-4:end,end-4:end);
prob.D = chol(prob.Q);

for i = 1:size(data,1)
    [sigma2(i,:),~,~,~] =kde(data(i,:),size(data,2),...
        min(data(i,:)),max(data(i,:)));
end

disp('Computing Moments for Cost')
tic
% E[w] -- mean
prob.muvec =  diracMixtureCostmean(data);
% E[W'W]
m2 = diracMixtureCostcov(data,sigma2);
prob.muvec2 = m2- prob.muvec.^2;
toc

disp('Solving ECF-STOC')
tic
[ECFSTOC_time_to_solve,ECFSTOC_total_time,ECFSTOC_opt_input_vector,...
    ECFSTOC_opt_mean_X,ECFSTOC_opt_val] = ECFSTOC(prob);
toc

for run_indx = 1:3
    fprintf('Run : %d\n',run_indx);
    [blackmore_time_to_solve(:,run_indx),blackmore_total_time(:,run_indx),blackmore_opt_input_vector(:,run_indx),...
    blackmore_opt_mean_X(:,run_indx),blackmore_opt_val(:,run_indx)] = BlackmoreTRo11PC(prob);
end

blackmore_mean_time = sum(blackmore_time_to_solve)/length(blackmore_time_to_solve);


%% Monte-Carlo Validation
        n_mcarlo_sims = 1e5;  
        Q = diag(prob.Q);
        Q = repmat(Q,1,n_mcarlo_sims);
        xtarget_mcarlo = repmat(prob.Xd, 1, n_mcarlo_sims);
        collection_of_input_vectors = [ECFSTOC_opt_input_vector,blackmore_opt_input_vector];
        collection_of_costs = [ECFSTOC_opt_val,blackmore_opt_val];
        % Have a collection of string
        collection_of_method_strings = {'ECFSTOC',...
                                        'BlackmoreTRO11_A',...
                                        'BlackmoreTRO11_B',...
                                        'BlackmoreTRO11_C'};
    
    
% SReachTools for Monte-Carlo simulation
max_rel_abserror = 0.1;
fprintf('Desired P{Hx<=g}: %1.2f | Desired relative abserror in cost: %1.2f\n',1-prob.Delta,max_rel_abserror);
for input_vec_indx = 1:size(collection_of_input_vectors,2)
    U_vector = collection_of_input_vectors(:,input_vec_indx);
    method_str = collection_of_method_strings{input_vec_indx};
    true_cost = collection_of_costs(input_vec_indx);
    % This function returns the concatenated state vector stacked columnwise
    n = n_mcarlo_sims;
    for i = 1:size(prob.Gd(:,1:2:end),2)
%         Wvec(2*i-1:2*i,:) = [normrnd(0,1,n,1)';normrnd(0,0.075,n,1)';];
    %     data(i,:) = [exprnd(1,n,1);];
%         Wvec(2*i-1:2*i,:) = [exprnd(2,n,1)'; normrnd(1,0.05,n,1)'];
%         Wvec(2*i-1:2*i,:) = [gamrnd(0.5,1,[n,1])'; unifrnd(-0.05,0.05,n,1)'];
        Wvec(2*i-1:2*i,:) = [2*wblrnd(5,4,[n,1])'; gamrnd(5,1,n,1)'];
    %     data(i,:) = [gamrnd(9,0.5,[n,1]); unifrnd(0,10,n,1);];
    %     data(i,:) = [gamrnd(8,0.5,[n,1]); exprnd(1,n,1);];
    %     data = repmat(data,size(Gd,2));
    end
     X_mcarlo = prob.Ad*prob.x0+prob.Bd*U_vector+prob.Gd*Wvec;
%      X_mcarlo = X.getRealizations(n_mcarlo_sims);
    % all does it column-wise
    particlewise_result = all(prob.pbig*X_mcarlo <= prob.qbig);
    prob_estim(input_vec_indx) = sum(particlewise_result)/n_mcarlo_sims;
    cost_estim(input_vec_indx) = 1/n_mcarlo_sims*sum(sum((X_mcarlo-xtarget_mcarlo).^2.*Q))+U_vector'*prob.R*U_vector;
    relative_abserror_in_cost(input_vec_indx) = abs(cost_estim(input_vec_indx) - true_cost)/true_cost;
    if prob_estim(input_vec_indx) >= 1 - prob.Delta && relative_abserror_in_cost(input_vec_indx) <= max_rel_abserror
        fprintf('PASSD: %s : Monte-Carlo via %1.0e particles | P{Hx<=g} = %1.3f | RelErr Cost = %1.3f\n',...
                method_str,... 
                n_mcarlo_sims,...
                prob_estim(input_vec_indx),...
                relative_abserror_in_cost(input_vec_indx));
    else
        fprintf('ERROR: %s : Monte-Carlo via %1.0e particles | P{Hx<=g} = %1.3f | RelErr Cost = %1.3f\n',...
                method_str,... 
                n_mcarlo_sims,...
                prob_estim(input_vec_indx),...
                relative_abserror_in_cost(input_vec_indx));
    end
end



%% Plotting

Fig4 = figure('Units', 'points', ...
                'Position', [0, 0, width, 200]); 

subplot(2,1,1)
hold on
h1 = yline(-q(1),'r','LineWidth',plot_linewidth);
yline(q(2),'r','LineWidth',plot_linewidth)
h11 = plot(1,prob.x0(1),'.b','MarkerSize',plot_markersize);
h2 = plot(2:(prob.T+1),prob.Xd(1:5:end),'go','MarkerSize',...
    plot_markersize,'LineWidth',2);
h3 = plot(2:(prob.T+1),ECFSTOC_opt_mean_X(1:5:end),'md',...
    'LineWidth',1,'MarkerSize',plot_markersize);
xl = prob.Ad*prob.x0+prob.Bd*ECFSTOC_opt_input_vector+GdWl';
xu = prob.Ad*prob.x0+prob.Bd*ECFSTOC_opt_input_vector+GdWu';
xm = prob.Ad*prob.x0+prob.Bd*ECFSTOC_opt_input_vector+prob.Gd*Wvec;
xlp = plot(2:(prob.T+1),xl(1:5:end),'r','LineWidth',1);
xup = plot(2:(prob.T+1),xu(1:5:end),'r','LineWidth',1);
h4 = plot(2:(prob.T+1),blackmore_opt_mean_X(1:5:end,1),...
    'ks','MarkerSize',plot_markersize);

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');

% xlabel('Time Step, k')
ylabel('Altitude, ft')
% legend([h1 h11 h2 h3 h4],{'Target Tube',...
%     'Initial state','Target Trajectory',...
%     'ECF Stochastic Optimal Control',...
%     sprintf('Particle control (PC), %i Particles',prob.N)},...
%     'Location','southoutside','FontSize',plot_fontSize);
box on;
set(gca,'FontSize',plot_fontSize);
set(gca,'xtick',[]);
axis([1 10 85000 85200])

subplot(2,1,2)
hold on
h1 = yline(-q(3),'r','LineWidth',plot_linewidth);
yline(q(4),'r','LineWidth',plot_linewidth)
h11 = plot(1,prob.x0(2),'.b','MarkerSize',plot_markersize);
h2 = plot(2:(prob.T+1),prob.Xd(2:5:end),'go','MarkerSize',...
    plot_markersize,'LineWidth',2);
h3 = plot(2:(prob.T+1),ECFSTOC_opt_mean_X(2:5:end),'md',...
    'LineWidth',1,'MarkerSize',plot_markersize);
xl = prob.Ad*prob.x0+prob.Bd*ECFSTOC_opt_input_vector+GdWl';
xu = prob.Ad*prob.x0+prob.Bd*ECFSTOC_opt_input_vector+GdWu';
xm = prob.Ad*prob.x0+prob.Bd*ECFSTOC_opt_input_vector+prob.Gd*Wvec;
xlp = plot(2:(prob.T+1),xl(2:5:end),'r','LineWidth',1);
xup = plot(2:(prob.T+1),xu(2:5:end),'r','LineWidth',1);
h4 = plot(2:(prob.T+1),blackmore_opt_mean_X(2:5:end,1),...
    'ks','MarkerSize',plot_markersize);

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');

xlabel('Time Step, k')
ylabel('Velocity, ft/s')
axis([1 10 7650 7750])
legend([h1 h11 h2 h3 xlp h4],{'Target Tube',...
    'Initial state','Target Trajectory',...
    'ECF Stochastic Optimal Control',...
    'Confidence interval',...
    sprintf('Particle control (PC), %i Particles',prob.N)},...
    'Location','southoutside','FontSize',plot_fontSize);
box on;
set(gca,'FontSize',plot_fontSize);

function m = diracMixtureCostmean(data)

	t = 0;

    m = (1i)^(-1)*1/(size(data,2))*sum(1i*data,2);


end

function m2 = diracMixtureCostcov(data,sigma)

	t = 0;

    m2 = (1i)^(-2)*1/(size(data,2))*...
        sum(-sigma.^2+(1i*data).^2,2);

end


function cf = diracMixtureICC(t,data,sigma)

	t = reshape(t,length(t),1);
    cf_int = sum(1/size(data,2)*exp(1i *  t * data),2).*exp(-sigma*(t).^2/2);
    cf = cf_int; 

end