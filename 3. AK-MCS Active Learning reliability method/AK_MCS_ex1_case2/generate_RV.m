function S = generate_RV(probdata,analysisopt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATION OF RANDOM NUMBERS
% Exctract model data
marg = probdata.marg;
R = probdata.correlation;
parameter = probdata.parameter;
% Find number of random variables
nrv = size(probdata.correlation,1);
% Modify correlation matrix and perform Cholesky decomposition
Ro = mod_corr( probdata, R );
Lo = (chol(Ro))';
% Generate uncorrelated standard normal variables u
u = mvnrnd(zeros(size(marg,1),1),eye(nrv),analysisopt.Nsamples )';  
z = Lo * u;
lambda = parameter(:,3);
xi = parameter(:,4);
x = [];
for i=1:nrv
    %x = [x; exp(z(i,:)* xi(i) + lambda(i))]; % for lognormal
    x = [x; z(i,:) * parameter(i,2) + parameter(i,1)];
end
S = x';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end