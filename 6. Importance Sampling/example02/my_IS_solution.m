clc,clear; format compact; format shortG;

addpath(genpath(strcat(pwd,'\ferumcore'  )) );
rng(19911109)

% Define Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params       = struct;

% Definition of PDF
probdata.marg(1,:) = [ 1  6  1  5 0 0 0 0 0]; % R
probdata.marg(2,:) = [ 1  3  1  3 0 0 0 0 0]; % S

% Define the correlation matrix
probdata.correlation   = eye(2);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MCS
u0        = zeros(nrv,1); % origin in standard normal space
stdv1     = 1;

% Initialize
i              = 0;
n_samples      = 10^3;

failed_samples = [];
safe_samples   = [];
while i < n_samples

    i = i+1;
    
    % draw samples
    ui = u0 + stdv1*eye(nrv)*randn(nrv,1);

    % evaluate the limit state function at the real space
    gi   = G_FUN( ui );

    Ii   = gi < 0;

    q(i) = Ii;

    %%%% for plots
    if Ii 
        failed_samples = [failed_samples, ui];
    else
        safe_samples   = [safe_samples, ui];
    end

end

Pf_mcs = (1/n_samples) * sum(q)


% plots
set(0,'defaultAxesFontSize',12)
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
set(groot, 'defaultLegendInterpreter', 'latex')

f1 = figure;
set(gcf,'units','inches','position',[1,1,6,2.6])
tiledlayout(1,2);

xx1 = linspace(-5,5,100);
xx2 = linspace(-5,5,100);

[X1, X2] = meshgrid(xx1, xx2);

for i=1:size(X1,1)
    for j=1:size(X1,2)
        xx = [X1(i,j);X2(i,j)];  
        gxx(i,j) = G_FUN( xx );
    end
end

nexttile;

scatter(failed_samples(1,:), failed_samples(2,:),8, [1 0 0], 'filled');
hold on;

scatter(safe_samples(1,:), safe_samples(2,:),8, [0.6 0.6 0.6], 'filled'); 
hold on;

contour(X1, X2, gxx, [0, 0], 'k-','LineWidth', 1.5); 
hold on;

text(0.98,0.02,['$P_{MCS} = ', num2str(Pf_mcs,'%.3f'), '$'], ...
    'Units','normalized', ...
    'Interpreter','latex', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','bottom');

xline(0,'k--'); yline(0,'k--');

xlabel('r'); ylabel('s');
title('Crude Monte Carlo Sim.')
box on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importance sampling

% Define parameters of importance sampling density h(u)
stdv1     = 1;
sigma     = stdv1^2 * eye(nrv);
invsigma  = inv(sigma);

% Initialize
i         = 0;
n_samples = 10^3;
q         = zeros(n_samples,1);  

failed_samples = [];
safe_samples   = [];
while i < n_samples

    i = i+1;
    
    % draw samples
    ui = ustar + stdv1*eye(nrv)*randn(nrv,1);

    % evaluate the limit state function at the real space
    gi   = G_FUN( ui );

    Ii   = gi < 0;

    % get pdf value in the original PDF f(u) <- standard normal
    f_ui = (2*pi)^(-nrv/2) * exp(-0.5*(ui'*ui));

    % get pdf value in importance PDF h(u) <- Normal(ustar, sigma^2*I)
    h_ui = 1/sqrt((2*pi)^nrv * det(sigma)) * ...
           exp(-0.5*(ui-ustar)'*invsigma*(ui-ustar));

    q(i) = Ii * (f_ui/h_ui);

    %%%% for plots
    if Ii 
        failed_samples = [failed_samples, ui];
    else
        safe_samples   = [safe_samples, ui];
    end

end

Pf_is = (1/n_samples) * sum(q)


nexttile;

scatter(failed_samples(1,:), failed_samples(2,:), 8, [1 0 0], 'filled'); 
hold on;

scatter(safe_samples(1,:), safe_samples(2,:), 8, [0.6 0.6 0.6], 'filled'); 
hold on;

contour(X1, X2, gxx, [0, 0], 'k-','LineWidth', 1.5); 
hold on;

scatter(ustar(1),ustar(2), 40,'green', 'filled', 'Marker','diamond'); 
hold on;

plot([0,ustar(1)],[0,ustar(2)], Color='green',LineWidth=1.0); 
hold on;

text(0.98,0.02,['$P_{IS} = ', num2str(Pf_is,'%.3f'), '$'], ...
    'Units','normalized', ...
    'Interpreter','latex', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','bottom');


xline(0,'k--'); yline(0,'k--');

xlim([-5,5]); ylim([-5,5]);
xlabel('r'); ylabel('s');
title('Importance sampling')
box on;

exportgraphics(gcf,'figure.pdf','ContentType','vector');
exportgraphics(gcf,'figure.png','Resolution',300)

rmpath(genpath(strcat(pwd,'\ferumcore'))   );