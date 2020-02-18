function [ECFSTOC_time_to_solve,ECFSTOC_total_time,ECFSTOC_opt_input_vector,...
    ECFSTOC_opt_mean_X,ECFSTOC_opt_val] = ECFSTOC(prob)
%ECFSTOC Summary of this function goes here
%   Detailed explanation goes here


    Delta = prob.Delta;
    Xd = prob.Xd; 
    Q = prob.Q; 
    R = prob.R;
    D = prob.D;
    Ad = prob.Ad; 
    Bd = prob.Bd;
    Gd = prob.Gd;
    muvec = prob.muvec;
    muvec2 = prob.muvec2;
    
    ulim = prob.ulim;
    
    
    % Hyperplane constraints:
    
    pbig = prob.pbig; 
    qbig = prob.qbig;
    
    n_lin_const = size(pbig,1);
    
    % PWA underapproximation of logPhi: 
    
    pu_m = prob.pu_m;
    pu_c = prob.pu_c;
    xlb = prob.xlb;
    
    % Initial conditions: 
    
    x0 = prob.x0;





%% Optimization Problem: 
cvx_precision best
tstart = tic;
cvx_begin quiet

    variable U(size(Bd,2),1);
    variable d(n_lin_const, 1);
    
    minimize((Ad*x0+Bd*U-Xd+Gd*muvec)'*Q*(Ad*x0+Bd*U-Xd+Gd*muvec)+trace(Q*Gd*diag(muvec2)*Gd') + U'*R*U)
    
    subject to 
    
            for i = 1:n_lin_const
                
                for j = 1:length(pu_m{i})
                    
                    pu_m{i}(j)*(qbig(i)-pbig(i,:)*(Ad*x0+Bd*U)) +...
                                            pu_c{i}(j) >= 1 - d(i,1);
                    qbig(i)-pbig(i,:)*(Ad*x0+Bd*U) >= xlb(i);
                    
                    
                end
                
            end
            
            
            sum(d) <= Delta;
            d >= 0;
            abs(U) <= ulim;


t1 = toc(tstart);
cvx_end
t2 = toc(tstart);


        ECFSTOC_time_to_solve = t2 - t1;
        ECFSTOC_total_time = cvx_cputime;

        if strcmpi(cvx_status,'Solved') || strcmpi(cvx_status,'Inaccurate/Solved')
            ECFSTOC_opt_mean_X = Ad*x0+Bd*U+Gd*muvec;
            ECFSTOC_opt_val = cvx_optval;
            ECFSTOC_opt_input_vector = U;            
        else
            ECFSTOC_opt_mean_X = nan(size(Ad,1),1);
            ECFSTOC_opt_val = nan;
            ECFSTOC_opt_input_vector = nan(size(Bd,2),1);         
        end
end

