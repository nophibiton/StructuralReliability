%
% Example 8.7 - Reliability analysis of ductile frame (Page 229)
%
% Der Kiureghian, A. (2022). Structural and System Reliability. 
% Cambridge University Press. https://doi.org/10.1017/9781108991889
%
clc,clear; format compact; format shortG;

addpath(genpath(strcat(pwd,'\ferumcore'  )) );
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

% Define transformation using Nataf model
Ro  = mod_corr(probdata, probdata.correlation );
Lo  = (chol(Ro))';    % Cholesky decomposition
u2x = @(u) u_to_x(u, probdata,Lo);

% Define limit state function in G(T(u))
G_FUN1 = @(u) g_fun1(u2x(u),params);
G_FUN2 = @(u) g_fun2(u2x(u),params);
G_FUN3 = @(u) g_fun3(u2x(u),params);
LSF    = {G_FUN1,G_FUN2,G_FUN3};

nrv    = size(probdata.marg,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importance sampling

% Get ustar from FORM
ncomp     = 3;
betas     = zeros(ncomp,1);
ustars    = zeros(nrv,ncomp);
for j=1:3

    [formresults]  = form(j,probdata,analysisopt,gfundata,[],[]);
    ustars(:,j)    = formresults.dsptu;
    betas(j)       = formresults.beta1;

end

% Define parameters of importance sampling density h(u)
% h(u) = sum( weight_i * N(ustar_i, sigma*I) )
% weight_i is proportional to Phi(-beta_i)
stdv1       = 1;
sigma       = stdv1^2 * eye(nrv);
invsigma    = inv(sigma);
Pf_form     = normcdf(-betas);
weights     = Pf_form / sum(Pf_form);
cum_weights = cumsum(weights);


% Initialize
i         = 0;
n_samples = 10^4;
q         = zeros(n_samples,1);  

failed_samples = [];
safe_samples   = [];
while i < n_samples

    i   = i+1;
    
    % select from which importance sampling density to sample
    idx = find(rand <= cum_weights,1,'first');

    % draw samples
    ui = ustars(:,idx) + stdv1*eye(nrv)*randn(nrv,1);

    % evaluate the limit state function at the real space
    gj   = zeros(ncomp,1);
    for j=1:ncomp
        gj(j)  = LSF{j}( ui );
    end
    Ii   = any(gj < 0); % series event

    % get pdf value in the original PDF f(ui) <- standard normal
    f_ui = (2*pi)^(-nrv/2) * exp(-0.5*(ui'*ui));

    % get pdf value in importance PDF 
    % h(ui) <- sum( weight_i * Normal(ustar_i, sigma^2*I) )
    h_ui = 0;
    for j=1:ncomp
        h_ui = h_ui + ( weights(j) )* 1/sqrt((2*pi)^nrv * det(sigma)) * ...
                     exp(-0.5*(ui-ustars(:,j))'*invsigma*(ui-ustars(:,j)));
    end

    q(i) = Ii * (f_ui/h_ui);

end

Pf_is = (1/n_samples) * sum(q)


rmpath(genpath(strcat(pwd,'\ferumcore'))   );