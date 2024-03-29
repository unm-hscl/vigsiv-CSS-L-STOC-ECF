function [blackmore_time_to_solve,blackmore_total_time,blackmore_opt_input_vector,...
    blackmore_opt_mean_X,blackmore_opt_val] = BlackmoreTRo11PC(prob)
    %% Blackmore TRo 2011 Code to stay in a feasible set. 
    % Coder: Vignesh Sivaramakrishnan
    
        N = prob.N; 
        T = prob.T;
        ulimu = prob.ulimu;
        uliml = prob.uliml;
        Delta = prob.Delta;
        x0 = prob.x0;
        W = prob.realizations;
        xtarget = prob.Xd;
        hbig = prob.pbig; 
        gbig = prob.qbig;
        Ad = prob.Ad;
        Bd = prob.Bd; 
        Gd = prob.Gd;
        Q = prob.Q;
        R = prob.R; 
        if isfield(prob,'xterm')
            xterm = prob.xterm;
        end
        

    % System matrices: 

        disp('---------BlackmorePC11-----------')
        disp(' ')
        fprintf('No. of particles: %d\n',N);

        large_constant = 5000;

    % Vectorize the target trajectory for the optimization problem. 
        xtargetbig = repmat(xtarget,1,N);

    % Randomly generate the disturbance vector from the standard normal.
            
        
%     if T >= 40
%         
%         blackmore_opt_mean_X = nan(length(Gd*W),1);
%         blackmore_opt_val = nan;
%         blackmore_opt_input_vector = nan(size(Bd,2),1); 
%         blackmore_time_to_solve = nan;
%         blackmore_total_time = nan;
%         warning('NOTE: BlackmoreTRo11 takes a very long time for long time horizons!!')
%         return;
%     else


    %% Run optimization problem for an optimal control policy
    % We run an optimization problem to determine the control policy over the
    % time horizon T.
        Q = diag(Q);
        Q = repmat(Q,1,N);
        tstart = tic;
        cvx_clear
            cvx_precision BEST
        cvx_begin 
            variable U_vector(size(Bd,2),1);
            variable xBl(size(Ad,1),N);
            variable mean_X(size(Ad,1),1);
            variable d(N) binary;

            minimize (1/N*sum(sum((xBl-xtargetbig).^2.*Q))+U_vector'*R*U_vector);

            subject to
              mean_X == Ad*x0+ Bd*U_vector;

              xBl(1:end,1:N) == Gd*W+repmat(mean_X,1,N);

                U_vector <= ulimu;
                U_vector >= uliml;

              for i = 1:N
                  
                  hbig*xBl(:,i) - gbig <= large_constant*d(i);                  
              end
              1/N*sum(d) >= 1 - Delta;
              
            if isfield(prob,'xterm')
                mean_X(end-4:end,:) == xterm;
            end
              

        t1 = toc(tstart);
        cvx_end;
        t2 = toc(tstart);
        blackmore_time_to_solve = t2 - t1;
        blackmore_total_time = cvx_cputime;

        if strcmpi(cvx_status,'Solved')
            blackmore_opt_mean_X = mean(xBl,2);
            blackmore_opt_val = cvx_optval;
            blackmore_opt_input_vector = U_vector;            
        else
            blackmore_opt_mean_X = nan(length(mean_GdTimesW),1);
            blackmore_opt_val = nan;
            blackmore_opt_input_vector = nan(size(Bd,2),1);         
        end
%     end

        fprintf('Total CVX Run Time for %1i particles: %1.4f seconds\n',...
            N,cvx_cputime)
        disp('------------------------------------')
        fprintf('Total CVX Solve Time for %1i particles: %1.4f seconds\n'...
            ,N,blackmore_time_to_solve)

        d = full(d);
end
