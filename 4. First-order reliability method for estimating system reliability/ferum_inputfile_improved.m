tic
clc, clear; format shortG; format compact;
%
% Citation
%
% Zhou, W., Gong, C., & Hong, H. P. (2017). New Perspective on Application 
% of First-Order Reliability Method for Estimating System Reliability. 
% Journal of Engineering Mechanics, 143(9). 
% https://doi.org/10.1061/(ASCE)EM.1943-7889.0001280
%
% Example 1: System Reliability of Pressurized Pipelines
% Containing Multiple Corrosion Defects
%

rng(0);
addpath(genpath(strcat(pwd,'\ferumcore')) )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINITION OF RANDOM VARIABLES
SMTS = 517;
P0 = 6.0;
tn = 7.16;

mean_sigma_u1 = 1.09*SMTS;
mean_p1 = 1.05*P0;
mean_d1 = 0.25*tn;
mean_e1 = 1.10;
mean_sigma_u2 = 1.09*SMTS;
mean_p2 = 1.05*P0;
mean_d2 = 0.30*tn;
mean_e2 = 1.10;

probdata.marg(1,:) = [ 2   mean_sigma_u1  0.03*mean_sigma_u1  mean_sigma_u1 0 0 0 0 0];
probdata.marg(2,:) = [ 15  mean_p1        0.1*mean_p1         mean_p1 0 0 0 0 0]; 
probdata.marg(3,:) = [ 16  mean_d1        0.2*mean_d1         mean_d1 0 0 0 0 0]; 
probdata.marg(4,:) = [ 2   mean_e1        0.172*mean_e1       mean_e1 0 0 0 0 0]; 
probdata.marg(5,:) = [ 2   mean_sigma_u2  0.03*mean_sigma_u2  mean_sigma_u2 0 0 0 0 0];
probdata.marg(6,:) = [ 15  mean_p2        0.1*mean_p2         mean_p2 0 0 0 0 0]; 
probdata.marg(7,:) = [ 16  mean_d2        0.2*mean_d2         mean_d2 0 0 0 0 0]; 
probdata.marg(8,:) = [ 2   mean_e2        0.172*mean_e2       mean_e2 0 0 0 0 0]; 

% Definition of correlation matrix
probdata.correlation(1,:) = [1.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0];
probdata.correlation(2,:) = [0.0, 1.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0];
probdata.correlation(3,:) = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.5, 0.0];
probdata.correlation(4,:) = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.5];
probdata.correlation(5,:) = [0.3, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0];
probdata.correlation(6,:) = [0.0, 0.8, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
probdata.correlation(7,:) = [0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 1.0, 0.0];
probdata.correlation(8,:) = [0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 1.0];

% Determine the parameters,the mean and standard deviation associated with each random variable
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


% Other Data
femodel = 0;
randomfield.mesh = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rzz = mod_corr(probdata, probdata.correlation );
Lj = (chol(Rzz))';
iLj = inv(Lj);
nrv = 8;
zd = zeros(nrv,1);

index = 1:1:nrv;

alpha_vector = [];
beta = [];

% Define the index which are included in the limit state function
j = [1,2,3,4]; % index of included
temp = index; temp(j)=[];
cj = temp; % index of not included

% Define a new probdata based on the RVs that are included
probdata_j.marg = probdata.marg(j,:);
probdata_j.correlation = probdata.correlation(j,j);
probdata_j.parameter = probdata.parameter(j,:);

% Limit State function (gfun) data
gfundata(1).evaluator = 'basic';
gfundata(1).type = 'expression';
gfundata(1).parameter = 'no';
gfundata(1).expression = 'g1(x)';

% Run FORM analysis from FERUM for g1(x)
[formresults] = form(1,probdata_j,analysisopt,gfundata,femodel,randomfield);

beta1 = formresults.beta1
beta = [beta; formresults.beta1];

uj = formresults.dsptu

% convert design point uj into the original space
Rzjj = Rzz(j,j );
Rzcjj = Rzz(cj,j ); 

Ljj = (chol(Rzjj ) )';

zj = Ljj*uj;
zcj = Rzcjj*(inv(Rzjj) )*zj;

% reorder indices
zd(j) = zj; 

zd(cj) = zcj

ud = iLj*zd

alpha_vector = [alpha_vector ud/formresults.beta1];

% Define the index which are included in the limit state function
j = [5,6,7,8]; % index of included
temp = index; temp(j)=[];
cj = temp; % index of not included

% Define a new probdata based on the RVs that are included
probdata_j.marg = probdata.marg(j,:);
probdata_j.correlation = probdata.correlation(j,j);
probdata_j.parameter = probdata.parameter(j,:);

% Limit State function (gfun) data
gfundata(2).evaluator = 'basic';
gfundata(2).type = 'expression';
gfundata(2).parameter = 'no';
gfundata(2).expression = 'g2(x)';

% Run FORM analysis from FERUM for g2(x)
[formresults] = form(2,probdata_j,analysisopt,gfundata,femodel,randomfield);

beta2 = formresults.beta1
beta = [beta; formresults.beta1];

uj = formresults.dsptu

% convert design point uj into the original space
Rzjj = Rzz(j,j );
Rzcjj = Rzz(cj,j ); 

Ljj = (chol(Rzjj ) )';

zj = Ljj*uj;
zcj = Rzcjj*(inv(Rzjj) )*zj;

% reorder indices
zd(j) = zj; zd(cj) = zcj

ud = iLj*zd

% calculate system failure probability
alpha_vector = [alpha_vector ud/formresults.beta1];
cmatrix = alpha_vector'*alpha_vector


rmpath(genpath(strcat(pwd,'\ferumcore')) )
toc
