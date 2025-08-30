function [r,R_DS,res_check,iter_results]=gen_DS_solver(R,m)

n=length(R);
IC=rand(n,m);
lb=-ones(n,m); ub=ones(n,m);
options=optimset('Display','off','TolFun',10^-16,'LargeScale','off',...
    'MaxFunEvals',10^5,'MaxIter',10^4,'TolCon',10^-8,'GradConstr','on'); 

[r,fval,exitflag]=fmincon(@(r)residual_gen_DS(r,R),IC,[],[],[],[],...
    lb,ub,@(r)cons_gen_DS(r),options)

res_check=sum(r.^2,2);
if any(res_check>1);
    n_row = find(res_check>1);
    r_ng = r(n_row,:);
    num_pow = abs(fix(log10(sum(r_ng.^2,2)-1)))-1;
    r_new = fix(r_ng.*(10.^num_pow*ones(1,size(r_ng,2))))...
        ./(10.^num_pow*ones(1,size(r_ng,2)));
    r(n_row,:) = r_new;
    res_check=sum(r.^2,2);
end
R_DS=r*r';
R_DS=R_DS-diag(diag(R_DS))+eye(n);
norm_err=norm(R_DS-R);
iter_results=[norm_err fval exitflag];

end