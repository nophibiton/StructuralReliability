clc, clear; format shortG; format compact;

rng(0);

addpath(genpath(strcat(pwd,'\dace')));

%
% Echard, B., Gayton, N., & Lemaire, M. (2011). AK-MCS: An active learning 
% reliability method combining Kriging and Monte Carlo Simulation. 
% Structural Safety, 33(2), 145–154. 
% https://doi.org/10.1016/J.STRUSAFE.2011.01.002
%
% Example 1: Case 1 (k=6)


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

        [gmean, mse] = predictor(S, krgMdl); % MSE represents the variance
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
beta = -norminv(Pf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare data for plotting

xTrain(1:size(initial_xTrain,1),:) = [];

x1 = linspace(-8,8,200);
x2 = linspace(-8,8,200);
[X,Y] = meshgrid(x1,x2);

Z=zeros(size(X));
Z_bar=zeros(size(X));
for i=1:size(X,1)
    for j=1:size(X,2)
        Z(i,j) = g_func([X(i,j),Y(i,j)]);
        Z_bar(i,j) = predictor([X(i,j),Y(i,j)], krgMdl);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot samples used in the active learning

f1 = figure;
set(f1,'units','inches','position',[1,1,5,5]);
leg_title = {};

s = scatter(S(:,1),S(:,2),'go','filled'); hold on; alpha(s,0.2); leg_title{1} = 'All samples';
scatter(xTrain(:,1),xTrain(:,2),'bo','filled'); hold on; leg_title{2} = 'Added samples';
scatter(initial_xTrain(:,1),initial_xTrain(:,2),100*ones(1,size(initial_xTrain,1)), ...
    'red','filled','pentagram'); hold on; leg_title{3} = 'Initial DOE';

contour(X,Y,Z,[0,0],LineColor="black",LineWidth=1); leg_title{4} = 'Exact g(x)';
contour(X,Y,Z_bar,[0,0],LineColor="magenta",LineWidth=1,LineStyle='--'); leg_title{5} = 'Surrogate g(x)';

xlabel('$x_1$',Interpreter='latex');
ylabel('$x_2$',Interpreter='latex');
legend(leg_title,Location="southeast");
box on

% Save figure content to PDF file
exportgraphics(gcf, 'example1_k6.pdf');


rmpath(genpath(strcat(pwd,'\dace')));