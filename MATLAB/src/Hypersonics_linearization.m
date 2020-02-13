%% abyvinod-vigsiv-PWA_JCC: Quadcopter with Additive Disturbance
%
% REQUIRED DEPENDENCIES: - MATLAB Symbolic Toolbox


%% Housekeeping
% clc, clear, close all

%% Parameters: 

S = 1.700E1; 
rho_0 = 6.7429E-5;
h_0 = 8.5E4;
h_s = 2.1358E4; 

Calpha_L = 4.6773E0; 
Cdeltae_L = 7.6224E-1; 
C0_L = -1.8714E-2; 

Calpha2_D = 5.8224E0;
Calpha_D = -4.5315E-2; 
Cdelta_e2_D = 8.1993E-1; 
Cdelta_e_D = 2.7699E-4; 
C0_D = 1.0131E-2; 

z_T = 8.36E0; 
coverline = 1.7E1; 
Calpha2_Malpha = 6.2926E0;
Calpha_Malpha = 2.1335E0;
C0_Malpha = 1.8979E-1; 
c_e = -1.2897E0; 

beta1 = -3.7693E5;
beta2 = -3.7225E4; 
beta3 = 2.6814E4; 
beta4 = -1.7277E4; 
beta5 = 3.5542E4; 
beta6 = -2.4216E3; 
beta7 = 6.3785E3; 
beta8 = -1.0090E2;

m = 300;
g = 32.17405;
Iyy = 5;


syms  h V alpha theta Q Phi deltae

%% Drag variables: 

C_L = Calpha_L*alpha + C0_L;
C_D = Calpha2_D*alpha^2 + Calpha_D*alpha + C0_D;
C_Malpha = Calpha2_Malpha*alpha^2+Calpha_Malpha*alpha+C0_Malpha;
C_Mdeltae = c_e*deltae;

Calpha3_T = beta1*Phi + beta2; 
Calpha2_T = beta3*Phi + beta4; 
Calpha_T = beta5*Phi + beta6; 
C0_T = beta7*Phi + beta8;

%% Parameters in Dynamics: 

L = 0.5*rho_0*V^2*S*C_L; 
D = 0.5*rho_0*V^2*S*C_D; 
T = Calpha3_T*alpha^3+Calpha2_T*alpha^2+Calpha_T*alpha+C0_T;
M = z_T*T+0.5*rho_0*V^2*S*coverline*[C_Malpha+C_Mdeltae]; 

%% Dynamics: 

f(1) = V*sin(theta-alpha);
f(2) = 1/m*(T*cos(alpha)-D)-g*sin(theta-alpha);
f(3) = 1/(m*V)*(-T*sin(alpha)-L)+Q+g/V*cos(theta-alpha);
f(4) = Q;
f(5) = M/Iyy;

vars = {h,V,alpha,theta,Q,Phi,deltae};

for i = 1:length(f)
    
    for j = 1:length(vars)
    
        jc(i,j) = jacobian(f(i),vars{j});
        
    end
 
end

replacem = {85000,7702.0808,deg2rad(1.5153),deg2rad(1.5153),0,0.2514,deg2rad(11.4635)};

A = double(subs(jc(:,1:5),vars,replacem));
Ts = 0.25;
% A = Ts*A+eye(5);
B = double(subs(jc(:,6:7),vars,replacem));
% B = 0.25*B;



