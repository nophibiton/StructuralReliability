function x = trans_u2x(u,probdata, Lo)

% unpack
distparams = probdata.distparams;
Xdist      = probdata.Xdist;

nrv = numel(u);
z   = Lo * u;
x   = zeros(size(z));

for i=1:nrv
    switch lower(Xdist{i})
        case 'normal'
            mu    = distparams(i,1);
            sigma = distparams(i,2);

            x(i)  = z(i)*sigma + mu;

        case 'lognormal'
            xi     = distparams(i,4);
            lambda = distparams(i,3);

            x(i)   = exp(z(i)*xi + lambda) ;          

        case 'gumbel'
            nu     = distparams(i,3);  
            alpha  = distparams(i,4);
            
            % Xi = F^(-1)[normcdf(Zi)] 
            p      = normcdf( z(i) );
            x(i)   = -(1/alpha)*log( -log(p) ) + nu;

        case 'weibull'
            %eps    = distparams(i,3); % special case of T3S (eps=0)  
            u      = distparams(i,4);
            k      = distparams(i,5);

            % Xi = F^(-1)[normcdf(Zi)] 
            p      = normcdf( z(i) );
            x(i)   = u*( -log(1-p) )^(1/k);
    end
  
end


end