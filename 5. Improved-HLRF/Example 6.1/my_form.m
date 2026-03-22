% Performs reliability analysis using First-order reliability method
% (FORM). The algorithm used is the improved Hasofer-Lind Rackwitz
% Fiessler. This code solves Example 6.1 - Jointly Normal Random Variables 
% found in [1].
%
% Reference:
% [1] Der Kiureghian, A. (2022). Structural and System Reliability. 
% Cambridge University Press. https://doi.org/10.1017/9781108991889
%
clc,clear; format compact; format shortG;

% define limit state parameters
lsfparams       = struct;

% define limit state function
lsffunc         = @g_fun;

% define probability data
probdata       = struct;
probdata.Xdist = {    'normal',     'normal'};
probdata.Xmu   = [          10,           20];
probdata.Xstd  = [           2,            5];
probdata.R     = [1.0,0.5; ...
                  0.5,1.0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get distribution parameters
distparams          = get_dist_params(probdata);
probdata.distparams = distparams;

Ro   = modify_corr(probdata.R,probdata);
Lo   = (chol(Ro))';
iLo  = inv(Lo);

% set starting point
x    = probdata.Xmu(:);
u    = trans_x2u(x,probdata,iLo);  % transform x to u

% set iteration number
i    = 0;

% set tolerances
eps1     = 10^-3;
eps2     = 10^-3;
max_iter = 20;

hist = struct;
hist.x = [];
hist.u = [];

while true

    % transform u to x
    x   = trans_u2x(u,probdata, Lo);


    % compute jacobians at x = xi
    Jux = calc_jacobian(x,u,probdata,Lo,iLo);
    Jxu = inv(Jux);

    % compute limit state function and gradients
    isGrad   = true;
    [gx,dgx] = eval_lsf(x,isGrad,lsffunc,lsfparams,distparams);
    Gu       = gx;
    dGu      = (dgx'*Jxu)';

    % set scale parameter
    if i == 0
        Go = Gu;
    end

    % compute alpha
    alpha = -dGu/norm(dGu);

    hist.x = [hist.x,x];
    hist.u = [hist.u,u];

    % check convergence
    if abs(Gu/Go) < eps1 && norm(u- alpha'*u*alpha) < eps2
        exit_message = 'MPP found.';
        break;
    end
    if i == max_iter
        exit_message = 'Maximum iteration is reached.';
        break;
    end

    % determine search direction
    d = (Gu/norm(dGu) + alpha'*u)*alpha-u;

    % determine step size
    %lambda = 1.0;
    lambda = calc_step_size(u,Gu,dGu,d, lsffunc,lsfparams,distparams,probdata,Lo);

    % determine new iteration point
    unew = u + lambda*d;
    u    = unew;
    i    = i + 1;

end

alpha = -dGu/norm(dGu)
ustar = u
xstar = trans_u2x(ustar,probdata, Lo)
beta  = alpha'*ustar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
set(0,'defaultAxesFontSize',12)
% set(0,'defaultTextFontName','Times New Roman')
% set(0,'defaultAxesFontName','Times New Roman')
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
set(groot, 'defaultLegendInterpreter', 'latex')


f1 = figure;
set(gcf,'units','inches','position',[1,1,6,2.8])
tiledlayout(1,2);

ax1 = nexttile;

xx1 = linspace(4,12,100);
xx2 = linspace(0,50,200);

[X1, X2] = meshgrid(xx1, xx2);

for i=1:size(X1,1)
    for j=1:size(X1,2)
        xx = [X1(i,j),X2(i,j)];  
        isGrad   = false;
        gxx(i,j) = eval_lsf(xx,isGrad,lsffunc,lsfparams,distparams);
    end
end

contour(X1, X2, gxx, [0, 0], 'k-','LineWidth', 2); hold on;

plot(hist.x(1,1:end),hist.x(2,1:end),'k-',LineWidth=1.0, ...
    Marker='o',MarkerFaceColor='k', MarkerSize=2); hold on;

scatter(hist.x(1,end),hist.x(2,end),100,Marker="pentagram", ...
    MarkerEdgeColor='k',MarkerFaceColor='r'); hold on;

text(hist.x(1,1)-0.3, hist.x(2,1)+2.5, ...
    '$\textbf{x}_0=\textbf{M}_x$',Interpreter='latex'); hold on;

xlabel('$x_1$',Interpreter='latex');
ylabel('$x_2$',Interpreter='latex');


ax2 = nexttile;

uu1 = linspace(-3,2,100);
uu2 = linspace(-2,3,200);

[U1, U2] = meshgrid(uu1, uu2);

for i=1:size(U1,1)
    for j=1:size(U1,2)
        uu = [U1(i,j);U2(i,j)];
        xx = trans_u2x(uu,probdata, Lo);
        isGrad   = false;
        Guu(i,j) = eval_lsf(xx,isGrad,lsffunc,lsfparams,distparams);
    end
end

quiver(0, 0, alpha(1), alpha(2), 1.5*beta, ...
    'LineWidth', 1.5,'MaxHeadSize', 0.3); hold on;

contour(U1, U2, Guu, [0, 0], 'k-','LineWidth', 2); hold on;

plot(hist.u(1,1:end),hist.u(2,1:end),'k-',LineWidth=1.0, ...
    Marker='o',MarkerFaceColor='k', MarkerSize=2); hold on;

scatter(hist.u(1,end),hist.u(2,end),100,Marker="pentagram", ...
    MarkerEdgeColor='k',MarkerFaceColor='r'); hold on;

text(hist.u(1,1)+0.1, hist.u(2,1)+0.15, ...
    '$\textbf{u}_0$',Interpreter='latex'); hold on;

xline(0);
yline(0);

xlabel('$u_1$',Interpreter='latex');
ylabel('$u_2$',Interpreter='latex');