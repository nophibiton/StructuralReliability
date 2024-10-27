clc, clear; format compact; format shortG;
%
%
% Kang, W. H., & Song, J. (2010). Evaluation of multivariate normal 
% integrals for general systems by sequential compounding. Structural Safety, 
% 32(1), 35–41. https://doi.org/10.1016/J.STRUSAFE.2009.06.001
%
%
% Example 3.2 Parallel system consisting of 10 components 
%   with equal reliabilityindexes and equal correlation coefficients
%

m = -3:1:3;

scm_value = zeros(1,length(m));
exact_value = zeros(1,length(m));

for l=1:length(m)

    beta_i = m(l)*ones(1,10);
    rho = eye(10,10);
    for ii =1:10
        for jj=1:10
            if ii~=jj
            rho (ii,jj) = 1 - (abs(ii-jj)/(10-1)); 
            end
        end
    end

    
    Correlation_matrix = rho;
    beta_sys = scm(beta_i, Correlation_matrix,'parallel');
    scm_value(l) = normcdf(-beta_sys)

end

semilogy(m,scm_value,Marker="none",LineStyle="-",Color="blue")
hold on

legend('SCM',Location='southwest')
xlabel('Component Reliability Index, \beta')
ylabel('Probability of Parallel System P(E_{sys})')
set(gcf,'units','inches','position',[1,1,4,3]);
ylim([10^-10,10^0])

