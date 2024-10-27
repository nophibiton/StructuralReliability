function obj = obj_fun_intersect(x,p1,p2, p3)
%
% find correlation coefficient of intersection events
%

rho_1and2k = x;

beta_1 = p1(1);
beta_2 = p1(2);
beta_k = p1(3);
beta_1and2 = p1(4);

rho_12 = p2(1);
rho_1k = p2(2);
rho_2k = p2(3);

% using matlab buitlin functions to calculate multinormal probabilities
if p3 ==0
    R3 = [1,rho_12,rho_1k ; rho_12,1,rho_2k ; rho_1k,rho_2k,1];
    R2 = [1,rho_1and2k ; rho_1and2k,1];
    trivar_cdf = mvncdf([-beta_1,-beta_2,-beta_k],zeros(1,3),R3);
    bivar_cdf = mvncdf([-beta_1and2,-beta_k],zeros(1,2),R2);

    obj = abs(trivar_cdf - bivar_cdf);
elseif p3==1
    R2 = [1,rho_12 ; rho_12,1];
    R2_comp = [1,rho_1and2k ; rho_1and2k,1];
    
    bivar_cdf = mvncdf([-beta_1,-beta_2],zeros(1,2),R2);
    bivar_cdf_comp = mvncdf([-beta_1and2,-beta_k],zeros(1,2),R2_comp);

    obj = abs(bivar_cdf-bivar_cdf_comp);
end