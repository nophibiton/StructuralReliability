function distparams = get_dist_params(probdata)

% unpack
Xdist = probdata.Xdist;
Xmu   = probdata.Xmu;
Xstd  = probdata.Xstd;

nrv        = numel(Xmu);
distparams = zeros(nrv,5);

for i=1:nrv
    switch lower(Xdist{i})
        case 'normal'
            mu    = Xmu(i);
            sigma = Xstd(i);

            distparams(i,1) = mu;  
            distparams(i,2) = sigma;
            distparams(i,3) = 0;  
            distparams(i,4) = 0;

        case 'lognormal'
            cov    = Xstd(i)/Xmu(i);
            zeta   = sqrt(log(1 + cov^2));
            lambda = log(Xmu(i)) - 0.5*zeta^2; 

            distparams(i,1) = Xmu(i);  
            distparams(i,2) = Xstd(i);
            distparams(i,3) = lambda;  
            distparams(i,4) = zeta;

        case 'gumbel'
            mu    = Xmu(i);
            std   = Xstd(i);
            euler_const = 0.5772156649; % eulers constant

            alpha = pi/(sqrt(6)*std);
            nu    = mu - (euler_const/alpha);

            distparams(i,1) = Xmu(i);  
            distparams(i,2) = Xstd(i);
            distparams(i,3) = nu;  
            distparams(i,4) = alpha;

        case 'weibull'
            mu    = Xmu(i);
            std   = Xstd(i);
            eps   = 0; % for weibull = TypeIII smallest (epsilon=0)

            cov = std / mu;
            k0 = cov^(-1.086);
            lambda0 = mu / gamma(1 + 1/k0);

            f = @(v)[eps + (v(1)-eps)*gamma(1+1/v(2))-mu ...
                     (v(1)-eps)*sqrt(gamma(1+2/v(2))-(gamma(1+1/v(2)))^2)-std ];
            
            x0 = [lambda0; k0]; 
            sol = fsolve(f, x0, optimoptions('fsolve','Display','off'));
            u = sol(1);
            k = sol(2);

            distparams(i,1) = Xmu(i);  
            distparams(i,2) = Xstd(i);
            distparams(i,3) = eps;  
            distparams(i,4) = u;
            distparams(i,5) = k;
    end

end

end