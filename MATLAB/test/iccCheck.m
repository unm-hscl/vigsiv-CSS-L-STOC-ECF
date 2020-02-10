%% vigsiv-CSS-L-STOC-ECF: iccCheck
% This script is used to check the validity of the ICC evaluation and the
% recover a underapproximation to a specific point, a. It utilizes double
% integrator dynamics and box constraints for the states.
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
    
n = 5000;
for i = 1:size(Gd,2)
%       data(i,:) = [normrnd(0,0.075,n,1);]';
%     data(i,:) = [exprnd(1,n,1);normrnd(1,0.05,n,1)]';
%     data(i,:) = [unifrnd(0.04,0.05,n,1); normrnd(0,0.04,n,1)]';
%     data(i,:) = [gamrnd(9,0.5,[n,1]); unifrnd(0.04,10,n,1);];
%     data(i,:) = [gamrnd(9,0.5,[n,1]); exprnd(1,n,1);];
%     data = repmat(data,size(Gd,2));
[bb(i,:),~,~,~] =kde(data,n,min(data(i,:)),max(data(i,:)));
end





% Generate the characteristic function to invert. 

    % Bounds on the safe set: 

    p = [1 0; -1 0;];
    q = linspace(5,3, prob.T);

    % Generate bounds for Ono and PWA: 
    pbig = kron(eye(prob.T),p);
    qbig = kron(q,[1,1])';

    n_lin_const = size(pbig,1);
    
    sigma = kron(eye(n_lin_const,length(bb)),[1,1]')*bb;
    
tic    
for k = 1:size(pbig,1)
                
                transform = pbig(k,:)*Gd;
                [ssigma(k),~,~,~] =kde(transform*data,n,min(transform*data),max(transform*data));
                cf_func = @(t) diracMixtureICC(t,data,transform,ssigma(k));
                clear options
                options.isPlot = false;
                options.xN = 1000; 
%                 options.N = 1000;
                result{k} = cf2DistGP(cf_func,[],[],options);
                
                x{k} = fliplr(result{1,k}.x)';
                pdf{k} = fliplr(result{1,k}.pdf)';
                cdf{k} = fliplr(result{1,k}.cdf)';
                
%                 x{k} = result{1,k}.x';
%                 pdf{k} = result{1,k}.pdf';
%                 cdf{k} = result{1,k}.cdf';

                [pu_m{k},pu_c{k},res(k,:)] =...
                    piecewiseUnder(x{k},cdf{k},1E-3,20);
                
                xind = find(x{k}==res(k,1));
                
                y{k} = min(pu_m{k}.*x{k}(xind:end)+pu_c{k},[],2);
                    
                wtransf = pbig(k,:)*Gd*data;
                figure(1)
                title('CDF')
                hold on
                empcdf = histogram(wtransf,'Normalization','cdf');
                hold on
                plot(x{k},cdf{k},'-b','LineWidth',2)
                plot(x{k}(xind:end),y{k},'-r','Linewidth',2)
                figure(2); 
                plot(x{k}(xind:end),cdf{k}(xind:end)-y{k})
                hold on
                
end
toc


function cf = diracMixtureICC(t,data,transform,sigma)

	t = reshape(t,length(t),1);
    cf_int = sum(1/size(data,2)*exp(1i *  t * transform * data),2).*exp(-sigma*(t).^2/2);
    cf = cf_int; 

end