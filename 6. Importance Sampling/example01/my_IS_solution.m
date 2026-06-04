clc,clear; format compact; format shortG;

addpath(genpath(strcat(pwd,'\ferumcore'  )) );
rng(19911109)

% base units
kN  = 1;
m   = 1;
kPa = kN/m^2;
MPa = 1000*kPa;

% Define Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params       = struct;
params.theta = 2;
params.A     = 0.190 *m^2;
params.S1    = 0.030 *m^3;
params.S2    = 0.015 *m^3;

% Definition of PDF
probdata.marg(1,:) = [ 1  250*kN*m  0.3*250*kN*m  250*kN*m 0 0 0 0 0]; % m1
probdata.marg(2,:) = [ 1  125*kN*m  0.3*125*kN*m  125*kN*m 0 0 0 0 0]; % m2
probdata.marg(3,:) = [ 15 2500*kN   0.2*2500*kN   2500*kN  0 0 0 0 0]; % p
probdata.marg(4,:) = [ 16 40*MPa    0.1*40*MPa    40*MPa   0 0 0 0 0]; % y

% Define the correlation matrix
probdata.correlation   = [1.0,0.5,0.3,0.0; ...
                          0.5,1.0,0.3,0.0; ...
                          0.3,0.3,1.0,0.0; ...
                          0.0,0.0,0.0,1.0];

probdata.parameter     = distribution_parameter(probdata.marg);


% Define models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define limit state function (gfun) data
gfundata(1).evaluator  = 'g_fun';
gfundata(1).parameter  = 'yes';
gfundata(1).thetag     = params;

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

[formresults]    = form(1,probdata,analysisopt,gfundata,[],[]);
ustar            = formresults.dsptu;

% Define transformation using Nataf model
Ro  = mod_corr(probdata, probdata.correlation );
Lo  = (chol(Ro))';    % Cholesky decomposition
u2x = @(u) u_to_x(u, probdata,Lo);

% Define limit state function in G(T(u))
G_FUN = @(u) g_fun(u2x(u),params); 


nrv       = size(probdata.marg,1);

% Define parameters of importance sampling density h(u)
stdv1     = 1;
sigma     = stdv1^2 * eye(nrv);
invsigma  = inv(sigma);

% Initialize
i         = 0;
n_samples = 10^4;
q         = zeros(n_samples,1);  
while i < n_samples

    i = i+1;
    
    % draw samples from h(u) importance density function
    ui = ustar + stdv1*eye(nrv)*randn(nrv,1);

    % evaluate the limit state function at the real space
    gi = G_FUN( ui );

    Ii = gi < 0;

    % get pdf value in the original PDF f(u) <- standard normal
    f_ui = (2*pi)^(-nrv/2) * exp(-0.5*(ui'*ui));

    % get pdf value in importance PDF h(u) <- Normal(ustar, sigma^2*I)
    h_ui = 1/sqrt((2*pi)^nrv * det(sigma)) * ...
           exp(-0.5*(ui-ustar)'*invsigma*(ui-ustar));

    q(i) = Ii * (f_ui/h_ui);

end

Pf = (1/n_samples) * sum(q)

rmpath(genpath(strcat(pwd,'\ferumcore'))   );