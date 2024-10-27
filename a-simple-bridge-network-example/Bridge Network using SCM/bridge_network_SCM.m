clc,clear; format compact; format shortG;

addpath(genpath(strcat(pwd,'\ferumcore')))

seed  = 1;
rng(seed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define statistics of random variables
% RANDOM VARIABLES
mean_R_TB = 0.4; std_R_TB = 0.3*mean_R_TB; % capacity R1,R2 Truss bridge
mean_R_GB = 0.5; std_R_GB = 0.5*mean_R_GB; % capacity R3,R4 Girder bridge
mean_R_CB = 0.3; std_R_CB = 0.3*mean_R_CB; % capacity R5 Cable bridge
mean_Si = 1.0; std_Si = 0.2*mean_Si; % site effect factor, Si
mean_A = 0.15; std_A = 0.5*mean_A; % peak ground acceleration, A
% Definition of PDF
probdata.marg(1,:) = [ 2  mean_R_TB  std_R_TB  mean_R_TB 0 0 0 0 0]; % R1
probdata.marg(2,:) = [ 2  mean_R_TB  std_R_TB  mean_R_TB 0 0 0 0 0]; % R2
probdata.marg(3,:) = [ 2  mean_R_GB  std_R_GB  mean_R_GB 0 0 0 0 0]; % R3
probdata.marg(4,:) = [ 2  mean_R_GB  std_R_GB  mean_R_GB 0 0 0 0 0]; % R4
probdata.marg(5,:) = [ 2  mean_R_CB  std_R_CB  mean_R_CB 0 0 0 0 0]; % R5
probdata.marg(6,:) = [ 2  mean_Si  std_Si  mean_Si 0 0 0 0 0]; % S
probdata.marg(7,:) = [ 2  mean_A  std_A  mean_A 0 0 0 0 0]; % A
% Definition of correlation matrix
probdata.correlation = eye(7);
% R1 and R2 are correlated
probdata.correlation(1,2) = 0.3; probdata.correlation(2,1) = 0.3;
% R3 and R4 are correlated
probdata.correlation(3,4) = 0.5; probdata.correlation(4,3) = 0.5;
% other parameters
probdata.parameter = distribution_parameter(probdata.marg);
% Analysis Options
analysisopt.ig_max    = 100;
analysisopt.il_max    = 5;
analysisopt.e1        = 0.001;
analysisopt.e2        = 0.001; 
analysisopt.step_code = 0;
analysisopt.grad_flag = 'FFD';
analysisopt.sim_point = 'dspt';
analysisopt.stdv_sim  = 1;
analysisopt.num_sim   = 100000;
analysisopt.target_cov = 0.0125;
% Find number of random variables
nrv = size(probdata.correlation,1);
% Other Data
femodel = 0;
randomfield.mesh = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = [];
Pf_comp = [];
beta_comp = [];

% define properties of limit state function
gfundata(1).evaluator = 'basic';
gfundata(1).type = 'expression';
gfundata(1).parameter = 'no';

% Limit State function (gfun) data for bridge 1
gfundata(1).expression = 'x(1)-x(6)*x(7)';

% Run modified FORM analysis from FERUM
ferum_form;
beta_comp = [beta_comp, formresults.beta1];
Pf_comp = [Pf_comp, formresults.pf1];
alpha = [alpha, formresults.alpha];

% Limit State function (gfun) data for bridge 2
gfundata(1).expression = 'x(2)-x(6)*x(7)';

% Run modified FORM analysis from FERUM
ferum_form;
beta_comp = [beta_comp, formresults.beta1];
Pf_comp = [Pf_comp, formresults.pf1];
alpha = [alpha, formresults.alpha];

% Limit State function (gfun) data for bridge 3
gfundata(1).expression = 'x(3)-x(6)*x(7)';

% Run modified FORM analysis from FERUM
ferum_form;
beta_comp = [beta_comp, formresults.beta1];
Pf_comp = [Pf_comp, formresults.pf1];
alpha = [alpha, formresults.alpha];


% Limit State function (gfun) data for bridge 4
gfundata(1).expression = 'x(4)-x(6)*x(7)';

% Run modified FORM analysis from FERUM
ferum_form;
beta_comp = [beta_comp, formresults.beta1];
Pf_comp = [Pf_comp, formresults.pf1];
alpha = [alpha, formresults.alpha];

% Limit State function (gfun) data for bridge 5
gfundata(1).expression = 'x(5)-x(6)*x(7)';

% Run modified FORM analysis from FERUM
ferum_form;
beta_comp = [beta_comp, formresults.beta1];
Pf_comp = [Pf_comp, formresults.pf1]
alpha = [alpha, formresults.alpha];

Ncomp = 5;

% correlation coefficient matrix
R = alpha'*alpha;

% SCM solution

% system event definition is a link set
% (E1 U E2)(E3 U E4)(E5)
sys_def = {[0 1 2 0 3 4 0 5 0],'link'};

beta_sys = scm(beta_comp, R, 'general',sys_def);
Pf_sys = normcdf(-beta_sys)

rmpath( genpath(strcat(pwd,'\ferumcore')) )
