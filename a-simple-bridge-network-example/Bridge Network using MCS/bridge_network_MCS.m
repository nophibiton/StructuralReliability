clc,clear; format compact; format shortG;

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
probdata.marg(6,:) = [ 2  mean_Si  std_Si  mean_Si 0 0 0 0 0]; % Si
probdata.marg(7,:) = [ 2  mean_A  std_A  mean_A 0 0 0 0 0]; % A
% Definition of correlation matrix
probdata.correlation = eye(7);
% R1 and R2 are correlated
probdata.correlation(1,2) = 0.3;
probdata.correlation(2,1) = 0.3;
% R3 and R4 are correlated
probdata.correlation(3,4) = 0.5;
probdata.correlation(4,3) = 0.5;
% other parameters
probdata.parameter = distribution_parameter(probdata.marg);
analysisopt.stdv_sim  = 1;
analysisopt.target_cov = 0.05;
% Find number of random variables
nrv = size(probdata.correlation,1);
% Exctract model data
marg = probdata.marg;
R = probdata.correlation;
parameter = probdata.parameter;
% Modify correlation matrix and perform Cholesky decomposition
Ro = mod_corr( probdata, R );
Lo = (chol(Ro))';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% (1) generate random variables
Nsamples = 1000000;
analysisopt.Nsamples   = Nsamples;
u = mvnrnd(zeros(size(marg,1),1),eye(nrv),analysisopt.Nsamples )';  

% (2) transform u to x using Nataf trasformation
z = Lo * u;
x = zeros(nrv, analysisopt.Nsamples );
for i=1:nrv
    % lognormal random variable
    xi = parameter(i,4 );
    lambda = parameter(i,3 );
    x(i,:) = exp(  z(i,:) * xi + lambda  );

end

% (3) evaluate each sample
sys_failure = 0;
comp1 = 0;
comp2 = 0;
comp3 = 0;
comp4 = 0;
comp5 = 0;
for k=1:Nsamples
    Ri = x(1:5,k);
    Si = x(6,k);
    A = x(7,k);

    % limit state function g(X) = Ri - A*Si for i=1,.. number of briges
    gi = Ri-A*Si;
    Ei = gi < 0;

    % system event definition is a link set
    % (E1 U E2)(E3 U E4)(E5)
    Esys = (Ei(1) | Ei(2)) & (Ei(3) | Ei(4)) & Ei(5);
    if Esys, sys_failure = sys_failure+1; end
    k

    if Ei(1), comp1 = comp1+1; end
    if Ei(2), comp2 = comp2+1; end
    if Ei(3), comp3 = comp3+1; end
    if Ei(4), comp4 = comp4+1; end
    if Ei(5), comp5 = comp5+1; end

end
Pf_comp = [comp1,comp2,comp3,comp4,comp5]/Nsamples
Pf_sys = sys_failure/Nsamples
