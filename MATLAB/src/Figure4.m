%% vigsiv-CSS-L-STOC-ECF: Figure 4: Double Integrator
% This code runs a simple double-integrator, LTI system with UNKNOWN 
% additive disturbance to evaluate the probability that the state will lie
% in a specified halfspace using approximate cdfs from the 
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
plot_markersize = 8;
plot_fontSize = 9;
plot_linewidth = 2;

%% System dynamics: 

    dim = 2; 
    prob.Ts = 0.25;
    prob.T = 10; % time horizon
    ulim = 0;

    input_space = Polyhedron('lb',-ulim,'ub',-ulim);
        disturb = RandomVector('UserDefined',...
                        @(N) [zeros(2,N)]);
    sys = getChainOfIntegLtiSystem(2, prob.Ts, input_space, disturb);
    sys_lti_no_input = LtiSystem('StateMatrix', sys.state_mat,...
                'DisturbanceMatrix', sys.dist_mat,...
                'Disturbance', sys.dist);
            
    % Generate concat matrices for evaluating the reach probability by
    % characteristic functions. 

    [prob.Ad,...
        prob.Bd,...
            prob.Gd] = getConcatMats(sys, prob.T); % Generate concat matrices.
    
% Generate data vector: 
    
n = 1000;
for i = 1:size(prob.Gd(:,1:dim:end),2)
%     data(2*i-1:2*i,:) = [normrnd(0,1,n,1)';normrnd(0,0.075,n,1)';];
%       data(i,:) = [exprnd(1,n,1);];
%     data(2*i-1:2*i,:) = [exprnd(2,n,1)'; normrnd(1,0.05,n,1)'];
%     data(2*i-1:2*i,:) = [gamrnd(0.5,1,[n,1])'; unifrnd(-0.05,0.05,n,1)'];
    data(2*i-1:2*i,:) = [unifrnd(-5,5,n,1)';0.005*gamrnd(8,0.5,[n,1])'; ];
%     data(i,:) = [gamrnd(9,0.5,[n,1]); unifrnd(0,10,n,1);];
%     data(i,:) = [gamrnd(8,0.5,[n,1]); exprnd(1,n,1);];
%     data = repmat(data,size(Gd,2));
    for j = 2*i-1:2*i
    [sigma(j),~,~,~] = kde(data(j,:),n,min(data(j,:)),max(data(j,:)));
    end
end

prob.N = 50;
for i = 1:size(prob.Gd(:,1:dim:end),2)
%     data(2*i-1:2*i,:) = [normrnd(0,1,n,1)';normrnd(0,0.075,n,1)';];
%       data(i,:) = [exprnd(1,n,1);];
%     data(2*i-1:2*i,:) = [exprnd(2,n,1)'; normrnd(1,0.05,n,1)'];
%     prob.realizations(2*i-1:2*i,:) = [gamrnd(0.5,1,[prob.N,1])'; unifrnd(-0.05,0.05,prob.N,1)'];
    prob.realizations(2*i-1:2*i,:) = [unifrnd(-5,5,prob.N,1)';0.005*gamrnd(8,0.5,[prob.N,1])'; ];
%     data(i,:) = [gamrnd(9,0.5,[n,1]); unifrnd(0,10,n,1);];
%     data(i,:) = [gamrnd(8,0.5,[n,1]); exprnd(1,n,1);];
%     data = repmat(data,size(Gd,2));
end

% Compute upper and lower bounds on the mean: 

for i = 1:size(data,1)
    
   
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

% Generate the characteristic function to invert. 

    % Bounds on the safe set: 

    p = [1 0; -1 0;];
    q = linspace(50,30, prob.T);

    % Generate bounds for Ono and PWA: 
    prob.pbig = kron(eye(prob.T),p);
    prob.qbig = kron(q,[1,1])';

    n_lin_const = size(prob.pbig,1);
disp('Generating PWAUA for ICC')    
tic    
for k = 1:n_lin_const
                
                transformed_rv = prob.pbig(k,:)*prob.Gd*data;
%                 [sigma(k),~,~,~] =kde(transformed_rv,n,min(transformed_rv),max(transformed_rv));
                cf_func = @(t) diracMixtureICC(t,transformed_rv,(prob.pbig(k,:)*prob.Gd)*diag(sigma)*(prob.pbig(k,:)*prob.Gd)');
                clear options
                options.isPlot = false;
                options.xN = 1000;
                options.xMin = min(transformed_rv); 
                options.xMax = max(transformed_rv);
%                 options.N = 1000;
                result{k} = cf2DistGP(cf_func,[],[],options);
                
                x{k} = fliplr(result{1,k}.x)';
                cdf{k} = fliplr(result{1,k}.cdf)';
                
%                 x{k} = result{1,k}.x';
%                 pdf{k} = result{1,k}.pdf';
%                 cdf{k} = result{1,k}.cdf';

                [prob.pu_m{k},prob.pu_c{k},res(k,:)] =...
                    piecewiseUnder(x{k},cdf{k},1E-3,20);
                
                prob.pu_m{k} = [prob.pu_m{k} 0];
                prob.pu_c{k} = [prob.pu_c{k} cdf{k}(end)];
                
                xind = find(x{k}==res(k,1));
                
                y{k} = min(prob.pu_m{k}.*x{k}(xind:end)+prob.pu_c{k},[],2);
                    
%                 figure(2)
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
prob.x0 = [0; 0];
Xd = linspace(-50,-50,prob.T)';
prob.Xd = reshape([Xd'; zeros(size(Xd'))],[],1);
prob.ulimu = 100;
prob.uliml = -100;
prob.xlb = res(:,1);

%% Generate Moments for the Quadratic Cost: 
prob.R = kron(eye(prob.T),0.01);
prob.Q = 10*eye(size(prob.Gd,2));
prob.D = chol(prob.Q);

disp('Computing Moments for Cost')
tic
% E[w] -- mean
prob.muvec =  diracMixtureCostmean(data);
% E[W'W]
m2 = diracMixtureCostcov(data,sigma');
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
    for i = 1:size(prob.Gd(:,1:dim:end),2)
%         Wvec(2*i-1:2*i,:) = [normrnd(0,1,n,1)';normrnd(0,0.075,n,1)';];
    %     data(i,:) = [exprnd(1,n,1);];
%         Wvec(2*i-1:2*i,:) = [exprnd(2,n,1)'; normrnd(1,0.05,n,1)'];
%         Wvec(2*i-1:2*i,:) = [gamrnd(0.5,1,[n,1])'; unifrnd(-0.05,0.05,n,1)'];
        Wvec(2*i-1:2*i,:) = [unifrnd(-5,5,n,1)';0.005*gamrnd(8,0.5,[n,1])'; ];
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
                'Position', [0, 0, width, height]); 
polyvertex =[0:prob.T,0:prob.T;[-prob.qbig(1),-prob.qbig(1:2:end)'],...
        [prob.qbig(1),prob.qbig(1:2:end)']]'; % Note: This needs MPT to run!!
P = Polyhedron('V',polyvertex);
hold on
h1 = plot(P,'alpha',0.1);
h11 = plot(0,prob.x0(1),'.b','MarkerSize',plot_markersize);
h2 = plot(1:(prob.T),Xd(1:end),'go','MarkerSize',...
    plot_markersize,'LineWidth',2);
h3 = plot(1:(prob.T),ECFSTOC_opt_mean_X(1:2:end),'md',...
    'LineWidth',1,'MarkerSize',plot_markersize);
xl = prob.Ad*prob.x0+prob.Bd*ECFSTOC_opt_input_vector+GdWl';
xu = prob.Ad*prob.x0+prob.Bd*ECFSTOC_opt_input_vector+GdWu';
xm = prob.Ad*prob.x0+prob.Bd*ECFSTOC_opt_input_vector+prob.Gd*Wvec;
xlp = plot(1:(prob.T),xl(1:2:end),'r','LineWidth',1);
xup = plot(1:(prob.T),xu(1:2:end),'r','LineWidth',1);
% for h = 1:n_mcarlo_sims/100
%     xmp = plot(2:(prob.T+1),xm(1:2:end,h));
% end
h4 = plot(1:(prob.T),blackmore_opt_mean_X(1:2:end,1),...
    'ks','MarkerSize',plot_markersize);


set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');

xlabel('Time Step, k')
ylabel('Position, x')
legend([h1 h11 h2 h3 xlp h4],{'Safety boundary',...
    'Initial state','Target trajectory',...
    'Algorithm 2',...
    'Confidence interval',...
    sprintf('Particle control, %i Particles',prob.N)},...
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