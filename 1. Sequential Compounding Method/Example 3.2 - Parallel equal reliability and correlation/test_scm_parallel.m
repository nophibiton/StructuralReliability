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
    rho = 0.5 * ones(10,10);
    for i = 1:size(rho,1)
        rho(i,i) = 1;
    end
    
    Correlation_matrix = rho;
    beta_sys = scm(beta_i, Correlation_matrix,'parallel');
    scm_value(l) = normcdf(-beta_sys)
    
    % exact value
    fun = @(s,beta,N) ((normcdf((-beta-sqrt(0.5).*s)./(sqrt(1-0.5)))).^(N)) .* normpdf(s) ;
    exact_value(l) = integral(@(s) fun(s,beta_i(l), 10),-inf,inf)

end

semilogy(m,scm_value,Marker="none",LineStyle="-",Color="blue")
hold on
semilogy(m,exact_value,Marker="o",LineStyle="none",MarkerEdgeColor="red")
hold off
legend('SCM','Exact',Location='southwest')
xlabel('Component Reliability Index, \beta')
ylabel('Probability of Parallel System P(E_{sys})')
set(gcf,'units','inches','position',[1,1,4,3]);
ylim([10^-7,10^0])
