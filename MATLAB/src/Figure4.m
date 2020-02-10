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
height = 200;
plot_markersize = 15;
plot_fontSize = 8;
plot_linewidth = 2;

% 