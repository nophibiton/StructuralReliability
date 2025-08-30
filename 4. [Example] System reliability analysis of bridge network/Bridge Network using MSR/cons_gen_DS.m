function [c,ceq,gradc,gradceq]=cons_gen_DS(r)

c=sum(r.^2,2)-(1-10^-6);
ceq=[];
gradc = [];
for ii = 1:size(r,2);
    GC = 2*diag(r(:,ii));
    gradc = [gradc; GC];
end
gradceq=[];

end