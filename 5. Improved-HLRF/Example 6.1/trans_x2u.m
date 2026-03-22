function u = trans_x2u(x,probdata,iLo)

% unpack
distparams = probdata.distparams;
Xdist      = probdata.Xdist;

nrv = numel(x);
z   = zeros(nrv,1);

for i=1:nrv
    switch lower(Xdist{i})
        case 'normal'
            mu    = distparams(i,1);
            sigma = distparams(i,2);

            z(i)  = (x(i)-mu)/sigma;

        case 'lognormal'
            xi     = distparams(i,4);
            lambda = distparams(i,3);

            z(i)   = ( log(x(i)) - lambda ) / xi ;          

        case 'gumbel'
            nu     = distparams(i,3);  
            alpha  = distparams(i,4);

            z(i)   = norminv(exp( -exp(-alpha*(x(i)-nu))));

        case 'weibull'
            eps    = distparams(i,3);  
            u      = distparams(i,4);
            k      = distparams(i,5);

            z(i)   = norminv( 1 - exp( -((x(i)-eps)/(u-eps))^k) );
    end

end

u = iLo * z;

end


