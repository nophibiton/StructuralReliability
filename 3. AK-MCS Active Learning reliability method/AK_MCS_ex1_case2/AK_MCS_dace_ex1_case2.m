clc, clear; format shortG; format compact;

rng(0);

addpath(genpath(strcat(pwd,'\dace')));

%
% Echard, B., Gayton, N., & Lemaire, M. (2011). AK-MCS: An active learning 
% reliability method combining Kriging and Monte Carlo Simulation. 
% Structural Safety, 33(2), 145–154. 
% https://doi.org/10.1016/J.STRUSAFE.2011.01.002
%
% Example 1: Case 2 (k=7)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINITION OF RANDOM VARIABLES
mean_x1 = 0; % x1
mean_x2 = 0; % x2
% Definition of PDF
probdata.marg(1,:) = [ 1  mean_x1  1  mean_x1 0 0 0 0 0]; % normal x1, std=1
probdata.marg(2,:) = [ 1  mean_x2  1  mean_x1 0 0 0 0 0]; % normal x2, std=1
% Definition of correlation matrix
probdata.correlation(1,1:2) = [1.0 0.0];
probdata.correlation(2,1:2) = [0.0 1.0];
% Determine the parameters,the mean and standard deviation associated with the distribution of each random variable
probdata.parameter = distribution_parameter(probdata.marg);
% Define analysis options
analysisopt.Nsamples = 10^6;
analysisopt.target_cov = 0.05;
analysisopt.NsamplesBatch = 10^5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrv = size(probdata.marg,1); % number of random variables

% generate random samples
S = generate_RV(probdata,analysisopt);

% initialize design of experiments (DoE)
N1 = 12;%(nrv+1)*(nrv+2)/2;
initial_DOE = randperm(size(S,1),N1);
xTrain = S(initial_DOE,:);
initial_xTrain = xTrain;
S(initial_DOE,:) = [];

% evaluate limit state function to the initial DoE
Y = g_func(xTrain);

while true

    while true
        
        % create Kriging model
        theta = [25 25]; lob = [1e-1 1e-1]; upb = [100 100];
        [krgMdl, perf] = dacefit(xTrain, Y, @regpoly1, @corrgauss, theta, lob, upb);

        [gmean, mse] = predictor(S, krgMdl);
        gsd = sqrt(mse);
    
        % calculate learning function U
        lf = abs(gmean ./ gsd);
        [min_U,ind] = min(lf);
        
        % enrich design of experiments with the min U
        xTrain = [xTrain; S(ind,:)];
        Y = [Y; g_func(S(ind,:))];
    
        % remove the added sample point from original population S
        S(ind,:) = [];
    
        % check stopping criterion of learning
        if min_U>2.0; break; end
    
        disp(['No. of Function evaluation: ', num2str(size(xTrain,1)),' to create surrogate model.']);
    end
    
    % calculate the responses of the limit state function using the 
    % surrogate model over all samples
    [y_hat,~] = predictor([S;xTrain], krgMdl);
    
    % calculate the probability of failure
    n_MC = size([S;xTrain],1); % number of samples
    Pf = sum(y_hat<=0)/n_MC;
    
    % calculate estimated COV of Pf
    estimatedCOV = sqrt((1-Pf)/(Pf*n_MC));
    
    if estimatedCOV < analysisopt.target_cov, break; end

    % Enrich new samples when COV of MCS is not attained
    analysisopt.Nsamples = analysisopt.NsamplesBatch;
    newSamples = generate_RV(probdata,analysisopt);
    S = [S;newSamples];


end


Pf
estimatedCOV

rmpath(genpath(strcat(pwd,'\dace')));