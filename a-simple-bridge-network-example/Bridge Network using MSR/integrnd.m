function int_val=integrnd(beta,r,c_sys,sys_type,varargin)

for i=1:length(varargin{1,1});
    for j=1:length(varargin);
        if j==1; n=i; else n=1; end
        s(j)=varargin{j}(n);
    end
    s=s';
    p_s=normcdf((-beta-r*s)./sqrt(1-sum(r.^2,2)));
    switch lower(sys_type)
        case 'parallel'
            int_val(i)=prod(p_s)*mvnpdf(s);
        case 'series'
            int_val(i)=prod(1-p_s)*mvnpdf(s);
        case 'general'
            [p]=prob_vector(p_s);
            int_val(i)=c_sys'*p*mvnpdf(s);
    end
    clear s
end

end