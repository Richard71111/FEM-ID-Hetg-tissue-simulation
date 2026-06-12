function [parameters, initial_state] = Initial_LR1(total_area)
%INITIAL_LR1 Define LR1 constants and the resting initial state.
% Input: total cell membrane area total_area in um^2.
% Outputs: parameters parameter structure and initial state vector initial_state.

% total_area - total cell area in um^2

% Standard ionic concentrations for ventricular cells
parameters.K_o = 5.4;                  % mM
parameters.K_i = 145;                  % mM
parameters.Na_o = 140;                 % mM
parameters.Ca_o = 1.8;                 % mM
parameters.Na_i = 10;                  % mM
parameters.F = 96.5;                   % Faraday constant, coulombs/mmol
parameters.R = 8.314;                  % gas constant, J/K
parameters.T = 273+37;                 % absolute temperature, K 
parameters.RTF=(parameters.R*parameters.T/parameters.F); % mV
% parameters.sqrt=sqrt(parameters.K_o/5.4); 

parameters.PK = 1.66e-6;                 % permability of K 
parameters.PNa_K = 0.01833;              % permability ratio of Na to K

% factor of 1e-8 converts from mS/cm^2 to mS/um^2 to mS
parameters.Isi_max=0.09e-8*total_area;               % mS
parameters.IKp_max=0.0183e-8*total_area;
parameters.Ib_max=0.03921e-8*total_area;
parameters.INa_max=16e-8*total_area;
parameters.IK1_max=0.6047e-8*total_area;
parameters.IK_max=0.282e-8*total_area;


%the initial conditions
% % Initial Gate Conditions */

v_init=-84.5286 ;  % mV
m_init = 0.0017; % sodium current activation gate
h_init =   0.9832;  % sodium current fast inactivation
J_init = 0.995484;   %  slow inactivation
d_init = 0.000003; % Calcium activation gate
f_init =  1.0000 ;  % Calcium  inactivation gate
X_init =  0.0057 ; % activation gate
Ca_init = 0.0002; % Calcium, mM


initial_state=[ v_init m_init h_init J_init d_init f_init X_init Ca_init];
end
