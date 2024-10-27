clc, clear; format shortG; format compact;

rng(0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINITION OF RANDOM VARIABLES
mean_x1 = 0; % x1
mean_x2 = 0; % x2
% Definition of PDF
probdata.marg(1,:) = [ 1  mean_x1  1  mean_x1 0 0 0 0 0]; % normal x1, mean=0, std=0.3
probdata.marg(2,:) = [ 1  mean_x2  1  mean_x1 0 0 0 0 0]; % normal x2, mean=0, std=0.3
% Definition of correlation matrix
probdata.correlation(1,1:2) = [1.0 0.0];
probdata.correlation(2,1:2) = [0.0 1.0];
% Determine the parameters,the mean and standard deviation associated with the distribution of each random variable
probdata.parameter = distribution_parameter(probdata.marg);
%
initial_Nsamples = 10^6;
analysisopt.Nsamples = initial_Nsamples;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATION OF RANDOM NUMBERS
nrv = size(probdata.marg,1); % number of random variables
% generate random samples
S = generate_RV(probdata,analysisopt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monte Carlo Simulation
fail=0;
for i=1:analysisopt.Nsamples
    g = g_func(S(i,:));
    if g<0
        fail=fail+1;
    end
end
Pf = fail/size(S,1);

while true

    Nsamples = size(S,1);
    cov_Pf = sqrt((1-Pf)/((Nsamples-1)*Pf));
    if cov_Pf<= 0.02, break; end

    analysisopt.Nsamples = 1;
    
    S = [S; generate_RV(probdata,analysisopt)];
    g = g_func(S(end,:)) ;

    if g<0
        fail=fail+1;
    end
    
    Pf = fail/size(S,1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

