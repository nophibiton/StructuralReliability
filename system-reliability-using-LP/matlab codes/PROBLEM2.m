clc, clear; format compact;
set(0,'defaultAxesFontSize',12)
set(0,'defaultTextFontName','Times New Roman')
set(0,'defaultAxesFontName','Times New Roman')


L_i = 5:0.1:8;
Pf_sys_exact=zeros(1,length(L_i));

f1 = figure;
set(f1,'units','inches','position',[1,1,4.5,3]);
leg_title = {};

for k=1:length(L_i)

    n = 6;
    L = L_i(k);

    x = zeros(n,1);
    for i=1:length(x), x(i) = L/(n-i+1); end
    
    lambda = 0.01;
    beta = 10;
    
    F = @(x,lambda,beta) 1-exp(-lambda*x.^(beta));
    b_n = F(x,lambda,beta);
    
    BB = zeros(n,n);
    for i=1:n
        bn = b_n(i);
        temp = [];
        for j=1:n-(i-1)
            temp = [temp,(bn^j)/(factorial(j))];
        end
        if i~=1
            BB(i,i-1:end) = [1,temp];
        else
            BB(1,1:end) = temp;
        end
        clear bn temp
    end
    Pf_sys_exact(k) = factorial(n) * det(BB);


    % calculate LP bounds
    [C] = event_matrix(n);
    n_mece = 2^n;

    P_Ei = zeros(1,length(x));
    for i=1:length(x)
        P_Ei(i) = Fi(x(i),i);
    end

    EiEj = nchoosek(1:n,2);
    for m=1:size(EiEj,1)
        P_EiEj(m)= Fij(x,EiEj(m,1),EiEj(m,2));
    end


    c_sys = [1;zeros(n_mece-1,1)];

    % Solve the problem using CVX
    cvx_begin
        variable p(n_mece)      
        maximize(c_sys' * p)    
        subject to
            sum(p) == 1      
    
            for i=1:n_mece
                p(i) >= 0 
            end
    
            for i=1:n
                C(:,i)'*p == P_Ei(i)
            end       
            
    cvx_end

    LP_unicomp_UB(k) = cvx_optval;

    % Solve the problem using CVX
    cvx_begin
        variable p(n_mece)      % Define the optimization variable
        maximize(c_sys' * p)    % Define the objective function
        subject to
            sum(p) == 1      
    
            for i=1:n_mece
                p(i) >= 0 
            end
    
            for i=1:n
                C(:,i)'*p == P_Ei(i)
            end

            for i=1:length(P_EiEj)
                (C(:,EiEj(i,1)).*C(:,EiEj(i,2)))'*p == P_EiEj(i)
            end            
            
    cvx_end

    LP_bicomp_UB(k) = cvx_optval;

end

semilogy(L_i,Pf_sys_exact,'k-',LineWidth=2.5); hold on; leg_title{1} = 'Exact';
semilogy(L_i,LP_unicomp_UB,'k-.'); hold on; leg_title{2} = 'Unicomponent (cvx)';
semilogy(L_i,LP_bicomp_UB,'k--'); hold on; leg_title{3} = 'Bicomponent (cvx)';

L_i = 5.7:0.1:8;

for k=1:length(L_i)

    n = 6;
    L = L_i(k);

    x = zeros(n,1);
    for i=1:length(x), x(i) = L/(n-i+1); end

    % calculate LP bounds
    [C] = event_matrix(n);
    n_mece = 2^n;

    P_Ei = zeros(1,length(x));
    for i=1:length(x)
        P_Ei(i) = Fi(x(i),i);
    end

    EiEj = nchoosek(1:n,2);
    for m=1:size(EiEj,1)
        P_EiEj(m)= Fij(x,EiEj(m,1),EiEj(m,2));
    end

    c_sys = [1;zeros(n_mece-1,1)];

    % Solve the problem using CVX
    cvx_begin
        variable p(n_mece)      % Define the optimization variable
        minimize(c_sys' * p)    % Define the objective function
        subject to
            sum(p) == 1      
    
            for i=1:n_mece
                p(i) >= 0 
            end
    
            for i=1:n
                C(:,i)'*p == P_Ei(i)
            end

            for i=1:length(P_EiEj)
                (C(:,EiEj(i,1)).*C(:,EiEj(i,2)))'*p == P_EiEj(i)
            end            
            
    cvx_end

    LP_bicomp_LB(k) = cvx_optval;

end

semilogy(L_i,LP_bicomp_LB,'k--'); hold on; leg_title{4} = '';

L_i = 7.3:0.1:8;

for k=1:length(L_i)

    n = 6;
    L = L_i(k);

    x = zeros(n,1);
    for i=1:length(x), x(i) = L/(n-i+1); end

    % calculate LP bounds
    [C] = event_matrix(n);
    n_mece = 2^n;

    P_Ei = zeros(1,length(x));
    for i=1:length(x)
        P_Ei(i) = Fi(x(i),i);
    end

    EiEj = nchoosek(1:n,2);
    for m=1:size(EiEj,1)
        P_EiEj(m)= Fij(x,EiEj(m,1),EiEj(m,2));
    end

    c_sys = [1;zeros(n_mece-1,1)];

    % Solve the problem using CVX
    cvx_begin
        variable p(n_mece)      
        minimize(c_sys' * p)    
        subject to
            sum(p) == 1      
    
            for i=1:n_mece
                p(i) >= 0 
            end
    
            for i=1:n
                C(:,i)'*p == P_Ei(i)
            end      
            
    cvx_end

    LP_unicomp_LB(k) = cvx_optval;

end

semilogy(L_i,LP_unicomp_LB,'k-.'); hold off; leg_title{5} = '';


xlabel('Load, L')
ylabel('Failure Probability, P_f')

legend(leg_title,Location="southeast");
filename = fullfile(strcat(pwd,'\Plots\','daniels_result.png'));
exportgraphics(gcf,filename,'Resolution',2000);
