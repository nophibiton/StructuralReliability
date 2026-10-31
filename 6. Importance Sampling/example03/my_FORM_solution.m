%
% Example 8.7 - Reliability analysis of ductile frame (Page 229)
%
% Der Kiureghian, A. (2022). Structural and System Reliability. 
% Cambridge University Press. https://doi.org/10.1017/9781108991889
%
clc,clear; format compact; format shortG;

addpath(genpath(strcat(pwd,'\ferumcore')) );
addpath(genpath(strcat(pwd,'\genz')) );
rng(19911109)

% Define Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params       = struct;

% Definition of PDF
probdata.marg(1,:) = [ 2  150  150*0.2  150 0 0 0 0 0]; % M1 [kN-m], LN
probdata.marg(2,:) = [ 2  150  150*0.2  150 0 0 0 0 0]; % M2 [kN-m], LN
probdata.marg(3,:) = [ 2  150  150*0.2  150 0 0 0 0 0]; % M3 [kN-m], LN
probdata.marg(4,:) = [ 2  150  150*0.2  150 0 0 0 0 0]; % M4 [kN-m], LN
probdata.marg(5,:) = [ 2  150  150*0.2  150 0 0 0 0 0]; % M5 [kN-m], LN
probdata.marg(6,:) = [15   50   50*0.4   50 0 0 0 0 0]; % H [kN], Gumbel
probdata.marg(7,:) = [ 3   60   60*0.2  150 0 0 0 0 0]; % V [kN], Gamma

% Define the correlation matrix
probdata.correlation   = eye(7);

for ik=1:5
    for jk=(ik+1):5
        probdata.correlation(ik,jk)=0.3;
    end
end
probdata.correlation = max(probdata.correlation,probdata.correlation');
probdata.parameter   = distribution_parameter(probdata.marg);


% Define models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define limit state function (gfun) data
gfundata(1).evaluator  = 'basic';
gfundata(1).type       = 'matlabfile';
gfundata(1).func       = @g_fun1;
gfundata(1).parameter  = 'yes';
gfundata(1).thetag     = params;

gfundata(2).evaluator  = 'basic';
gfundata(2).type       = 'matlabfile';
gfundata(2).func       = @g_fun2;
gfundata(2).parameter  = 'yes';
gfundata(2).thetag     = params;

gfundata(3).evaluator  = 'basic';
gfundata(3).type       = 'matlabfile';
gfundata(3).func       = @g_fun3;
gfundata(3).parameter  = 'yes';
gfundata(3).thetag     = params;

% Define analysis options
analysisopt.ig_max     = 50;
analysisopt.il_max     = 5;
analysisopt.e1         = 0.001;
analysisopt.e2         = 0.001; 
analysisopt.step_code  = 0;
analysisopt.grad_flag  = 'ffd';
analysisopt.sim_point  = 'dspt';
analysisopt.stdv_sim   = 1;
analysisopt.num_sim    = 100000;
analysisopt.target_cov = 0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORM
ncomp     = 3;
nrv       = size(probdata.marg,1);
alpha_mat = zeros(nrv,ncomp);
beta_i    = zeros(ncomp,1);
for i=1:3
    [formresults]  = form(i,probdata,analysisopt,gfundata,[],[]);
    ustar          = formresults.dsptu;
    
    beta_i(i)      = formresults.beta1;
    alpha_mat(:,i) = formresults.alpha;

end

Rs = alpha_mat'*alpha_mat;
Pf = 1 - qsimvnv(20000, Rs, -inf*ones(1,ncomp), beta_i)

rmpath(genpath(strcat(pwd,'\ferumcore')) );
rmpath(genpath(strcat(pwd,'\genz')) );