function obj = obj_fun_union(x,p1,p2,p3)
%
% find correlation coefficient of union events
%

rho_1or2k = x;

beta_1 = p1(1);
beta_2 = p1(2);
beta_k = p1(3);
beta_1or2 = p1(4);

rho_12 = p2(1);
rho_1k = p2(2);
rho_2k = p2(3);

if p3 ==0
    A = normpdf(-beta_k)/normcdf(-beta_k);
    B = A * (-beta_k + A);
    
    beta_1_k = (beta_1-rho_1k*A) / sqrt(1-B*(rho_1k)^2);
    beta_2_k = (beta_2-rho_2k*A) / sqrt(1-B*(rho_2k)^2);
    rho_12_k = (rho_12-rho_1k*rho_2k*B) / (sqrt(1-B*(rho_1k)^2)*sqrt(1-B*(rho_2k)^2));
    
    beta_1or2_k = (beta_1or2 - rho_1or2k*A) / sqrt(1-B*(rho_1or2k).^2);
    
    R2 = [1,rho_12_k;rho_12_k,1];
    obj = abs(1-mvncdf([beta_1_k,beta_2_k],zeros(1,2),R2)-normcdf(-beta_1or2_k));
elseif p3==1
    R2 = [1,rho_12 ; rho_12,1];
    R2_comp = [1,rho_1or2k ; rho_1or2k,1];
    
    bivar_cdf = mvncdf([beta_1,beta_2],zeros(1,2),R2);
    bivar_cdf_comp = mvncdf([rho_1or2k,beta_k],zeros(1,2),R2_comp);

    obj = abs(bivar_cdf-bivar_cdf_comp);
end

end