clc,clear; format compact; format shortG;

seed  = 1;
rng(seed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define statistics of random variables
% RANDOM VARIABLES
mean_R = 0.5; std_R = 0.5*mean_R; % capacity R 
mean_S = 1.0; std_S = 0.2*mean_S; % site effect factor, S
mean_A = 0.15; std_A = 0.5*mean_A; % peak ground acceleration, A
% Definition of PDF
probdata.marg(1,:) = [ 2  mean_R  std_R  mean_R 0 0 0 0 0]; % R1
probdata.marg(2,:) = [ 2  mean_S  std_S  mean_S 0 0 0 0 0]; % Si
probdata.marg(3,:) = [ 2  mean_A  std_A  mean_A 0 0 0 0 0]; % A
% Definition of correlation matrix
probdata.correlation = eye(3);
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
Nsamples = 10^7;
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
failure = 0;
for k=1:Nsamples
    Ri = x(1,k);
    Si = x(2,k);
    A = x(3,k);

    % limit state function g(X) = Ri - A*Si for i=1,.. number of briges
    gi = Ri-A*Si;
    Ei = gi < 0;

    if Ei, failure = failure + 1; end

end

Pf = failure/Nsamples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use normal random variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% transorm statistical properties from lognormal to normal
log_to_norm_var = @(log_cov) log(1+log_cov.^2);
log_to_norm_mean = @(log_mean,norm_var) log(log_mean)-0.5*norm_var;

log_mean_R = 0.5;
log_cov_R = 0.5;
norm_var_R = log_to_norm_var(log_cov_R);
norm_std_R = sqrt(norm_var_R);
norm_mean_R = log_to_norm_mean(log_mean_R,norm_var_R);

log_mean_S = 1.0;
log_cov_S = 0.2;
norm_var_S = log_to_norm_var(log_cov_S);
norm_std_S = sqrt(norm_var_S);
norm_mean_S = log_to_norm_mean(log_mean_S,norm_var_S);

log_mean_A = 0.15;
log_cov_A = 0.5;
norm_var_A = log_to_norm_var(log_cov_A);
norm_std_A = sqrt(norm_var_A);
norm_mean_A = log_to_norm_mean(log_mean_A,norm_var_A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define statistics of random variables
% RANDOM VARIABLES
% mean_R = 0.4; std_R = 0.3*mean_R; % capacity R 
% mean_S = 1.0; std_S = 0.2*mean_S; % site effect factor, S
% mean_A = 0.15; std_A = 0.5*mean_A; % peak ground acceleration, A
% Definition of PDF
probdata.marg(1,:) = [ 1  norm_mean_R  norm_std_R  norm_mean_R 0 0 0 0 0]; % R1
probdata.marg(2,:) = [ 1  norm_mean_S  norm_std_S  norm_mean_S 0 0 0 0 0]; % Si
probdata.marg(3,:) = [ 1  norm_mean_A  norm_std_A  norm_mean_A 0 0 0 0 0]; % A
% Definition of correlation matrix
probdata.correlation = eye(3);
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
% Nsamples = 10;
analysisopt.Nsamples   = Nsamples;
u = mvnrnd(zeros(size(marg,1),1),eye(nrv),analysisopt.Nsamples )';  

% (2) transform u to x using Nataf trasformation
z = Lo * u;
x = zeros(nrv, analysisopt.Nsamples );
for i=1:nrv
    % normal random variable
    x(i,:) = z(i,:)*parameter(i,2) + parameter(i,1 );

end

% (3) evaluate each sample
failure = 0;
Zi = [];
for k=1:Nsamples
    Ri = x(1,k);
    Si = x(2,k);
    A = x(3,k);

    % limit state function g(X) = Ri - A*Si for i=1,.. number of briges
    gi = Ri-A-Si;
    Ei = gi < 0;

    if Ei, failure = failure + 1; end

end

mean_Zi = mean(x(1,:) - x(2,:) - x(3,:) )
Pf = failure/Nsamples