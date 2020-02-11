%% vigsiv-CSS-L-STOC-ECF: iccCheck
% This script is used to check the validity of the moments derived for the
%  cost. It utilizes double integrator dynamics
%
% REQUIRED DEPENDENCIES: - CharFunTool 
%                          (https://github.com/witkovsky/CharFunTool/)
%                        - SReachTools
%                          (https://unm-hscl.github.io/SReachTools/)
%                        - MATLAB Statistics and Machine Learning
%                          Toolbox

%% Housekeeping
clc, clear, close all

%% System dynamics: 

    dim = 2; 
    prob.Ts = 0.25;
    prob.T = 10; % time horizon
    ulim = 0;
    uvec = 0;

    input_space = Polyhedron('lb',-ulim,'ub',-ulim);
        disturb = RandomVector('UserDefined',...
                        @(N) [zeros(2,N)]);
    sys = getChainOfIntegLtiSystem(2, prob.Ts, input_space, disturb);
    sys_lti_no_input = LtiSystem('StateMatrix', sys.state_mat,...
                'DisturbanceMatrix', sys.dist_mat,...
                'Disturbance', sys.dist);
            
    % Generate concat matrices for evaluating the reach probability by
    % characteristic functions. 

    [Ad, Bd, Gd] = getConcatMats(sys, prob.T); % Generate concat matrices.
    
% Generate data vector: 
    
n = 1000;
for i = 1:size(Gd,2)
      data(i,:) = [normrnd(1,0.075,n,1);]';
%       data(i,:) = [exprnd(1,n,1);];
%     data(i,:) = [exprnd(1,n,1); normrnd(1,0.05,n,1)]';
%     data(i,:) = [unifrnd(0.04,0.05,n,1); normrnd(0,0.04,n,1)]';
%     data(i,:) = [gamrnd(9,0.5,[n,1]); unifrnd(0.04,10,n,1);];
%     data(i,:) = [gamrnd(8,0.5,[n,1]); exprnd(1,n,1);];
%     data = repmat(data,size(Gd,2));
[sigma(i,:),~,~,~] =kde(data,n,min(data(i,:)),max(data(i,:)));
end

%% Generate Moments for the Quadratic Cost: 

Q = 1*eye(size(Gd,2));
D = chol(Q);

trans_data = D*Gd*data;

disp('Computing Moments')
tic
% D'*E[DGW] -- mean
m = D'*diracMixtureCostmean(data); 
% E[W'G'D'*DGW]
m2 = diracMixtureCostcov(data,sigma);
var = m2 - m.^2;
toc

function m = diracMixtureCostmean(data)

	t = 0;

    m = (1i)^(-1)*1/(size(data,2))*sum(1i*data,2);


end

function m2 = diracMixtureCostcov(data,sigma)

	t = 0;

    m2 = (1i)^(-2)*1/(size(data,2))*...
        sum(-sigma.^2+(1i*data).^2,2);

end