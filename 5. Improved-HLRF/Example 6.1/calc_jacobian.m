function Jux = calc_jacobian(x,u,probdata,Lo,iLo)

% unpack
distparams = probdata.distparams;
Xdist      = probdata.Xdist;

nrv = numel(u);
z   = Lo * u;

Jzx = zeros(nrv,nrv);

for i=1:nrv
    switch lower(Xdist{i})
        case 'normal'
            sigma = distparams(i,2);

            Jzx(i,i)  = 1/sigma;

        case 'lognormal'
            xi = sqrt( log( 1 + ( distparams(i,2) / distparams(i,1) )^2 ) );

            Jzx(i,i)  = 1/(xi*x(i));         
         
        case 'gumbel'
            nu     = distparams(i,3);  
            alpha  = distparams(i,4);
            
            % Jux = inv(Lo) * diag[ fi(xi)/normpdf(zi) ]
            pdf_xi    = alpha*exp(-alpha*(x(i)-nu)-exp(-alpha*(x(i)-nu)));
            Jzx(i,i)  = pdf_xi/normpdf( z(i) );

        case 'weibull'
            %eps    = distparams(i,3); % special case of T3S (eps=0)  
            u      = distparams(i,4);
            k      = distparams(i,5);

            % Jux = inv(Lo) * diag[ fi(xi)/normpdf(zi) ]
            pdf_xi    = (k/u)*(x(i)/u)^(k-1)*exp(-(x(i)/u)^k);
            Jzx(i,i)  = pdf_xi/normpdf( z(i) );
    end
  
end

Jux = iLo * Jzx;

end