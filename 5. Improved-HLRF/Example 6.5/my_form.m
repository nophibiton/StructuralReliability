% Performs reliability analysis using First-order reliability method
% (FORM). The algorithm used is the improved Hasofer-Lind Rackwitz
% Fiessler. This code solves Example 6.5 - Reliability of a column 
% found in [1].
%
% Reference:
% [1] Der Kiureghian, A. (2022). Structural and System Reliability. 
% Cambridge University Press. https://doi.org/10.1017/9781108991889
%
clc,clear; format compact; format shortG;

% define base units
kN  = 1;
m   = 1;
kPa = kN/m^2;
MPa = 1000*kPa;

% define limit state parameters
lsfparams       = struct;
lsfparams.theta = 2;
lsfparams.A     = 0.190 *m^2;
lsfparams.S1    = 0.030 *m^3;
lsfparams.S2    = 0.015 *m^3;

% define limit state function
lsffunc         = @g_fun;

% define probability data
probdata       = struct;
probdata.Xdist = {    'normal',     'normal',    'gumbel',  'weibull'};
probdata.Xmu   = [    250*kN*m,     125*kN*m,     2500*kN,     40*MPa];
probdata.Xstd  = [250*kN*m*0.3, 125*kN*m*0.3, 2500*kN*0.2, 40*MPa*0.1];
probdata.R     = [1.0,0.5,0.3,0.0; ...
                  0.5,1.0,0.3,0.0; ...
                  0.3,0.3,1.0,0.0; ...
                  0.0,0.0,0.0,1.0];

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

    % check convergence
    if abs(Gu/Go) < eps1 || norm(u- alpha'*u*alpha) < eps2
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
    % lambda = 1.0; % HL-RF algorithm
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
disp(exit_message)