clc, clear; format compact; format shortG;
%
%
% Kang, W. H., & Song, J. (2010). Evaluation of multivariate normal 
% integrals for general systems by sequential compounding. Structural Safety, 
% 32(1), 35–41. https://doi.org/10.1016/J.STRUSAFE.2009.06.001
%
%
% Example 3.6 Parallel system consisting of 5–50 components with equal
%    reliability indexes and equal correlation coefﬁcients

N = 5:5:50; 
scm_value = zeros(3,length(N));
exact_value = zeros(3,length(N));

for l=1:3
    for k=1:length(N)
    
        % using scm
        beta_i = l*ones(1,N(k));
        rho = 0.5 * ones(N(k),N(k));
        for i = 1:size(rho,1)
            rho(i,i) = 1;
        end
        Correlation_matrix = rho;
        beta_sys = scm(beta_i, Correlation_matrix,'parallel');
        scm_value(l,k) = normcdf(-beta_sys);
        
        % exact value
        fun = @(s,beta,N) ((normcdf((-beta-sqrt(0.5).*s)./(sqrt(1-0.5)))).^(N)) .* normpdf(s) ;
        exact_value(l,k) = integral(@(s) fun(s,beta_i(1),N(k)),-inf,inf);
    end

end

figure
for l=1:3
    semilogy(N,scm_value(l,:),Marker="none",LineStyle="-",Color="blue")
    hold on
    semilogy(N,exact_value(l,:),Marker="o",LineStyle="none",MarkerEdgeColor="red")
    text(N(6),scm_value(l,6)*2.5,strcat('\beta = ',num2str(l)))
end
hold off
legend('SCM','Exact')
ylim([10^-10 10^0])
xlim([5 50])
grid on
xlabel('No. of components, N')
ylabel('Probability of Parallel System P(E_{sys})')
set(gcf,'units','inches','position',[1,1,4,3]);
