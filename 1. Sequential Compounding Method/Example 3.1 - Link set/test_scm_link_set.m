clc, clear; format compact; format shortG;
%
%
% Kang, W. H., & Song, J. (2010). Evaluation of multivariate normal 
% integrals for general systems by sequential compounding. Structural Safety, 
% 32(1), 35–41. https://doi.org/10.1016/J.STRUSAFE.2009.06.001
%
%
% Example 3.1 Illustrative example: a link-set system with five components
%


beta_i = [-1;-1;-1;-1;-1]'

R = [1.0, 0.8, 0.6, 0.4, 0.2; ...
     0.8, 1.0, 0.8, 0.6, 0.4; ...
     0.6, 0.8, 1.0, 0.8, 0.6; ...
     0.4, 0.6, 0.8, 1.0, 0.8; ...
     0.2, 0.4, 0.6, 0.8, 1.0]

%
%  Ʌ : Intersection, V : Union
% 
%  Esys = E1 Ʌ E2 Ʌ (E3 V E4) Ʌ E5
sys_def = {[0 1 0 2 0 3 4 0 5 0],'link'};

beta_sys = scm(beta_i, R, 'general',sys_def);
Pfsys = normcdf(-beta_sys)
