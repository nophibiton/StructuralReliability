function mvnval  = mvnlps( mu, sigma, q, e, r, re )
%
% MVNLPS Multivariate Normal Distribution Value for an ellipsoid.        
%        MVNVAL = MVNLPS( MU, SIGMA, Q, E, R, RE ) computes
%        an MVN value to relative accuracy RE for an ellipsoid centered 
%        at Q with radius R and positive semi-definite ellipsoid matrix E: 
%                MVNVAL = PROB( ( X - Q )'E ( X - Q ) < R^2 )
%        SIGMA is a positive definite covariance matrix for a multivariate
%        normal (MVN) density with mean MU. MU and Q must be column vectors.
%        Example:
%         sg = [ 3 2 1;2 2 1;1 1 1]; mu = [1 -1 0]';        
%         e = [4 1 -1; 1 2 1; -1 1 2]; q = [2 3 -2]';
%         p = mvnlps( mu, sg, q, e, 4, 1e-5 ); disp(p)
%
%        Reference for basic algorithm is (S&O):
%         Sheil, J. and O'Muircheartaigh, I. (1977), Algorithm AS 106:
%           The Distribution of Non-negative Quadratic Forms in Normal
%           Variables, Applied Statistics 26, pp. 92-98.
%
%        Matlab implementation by
%         Alan Genz, WSU Math, Pullman, WA 99164-3113; alangenz@wsu.edu.
%
%  Transformation to problem with diagonal covariance matrix
%
   [ n, n ] = size(sigma); L = chol( sigma ); [ V D ]  = eig( L*e*L' );
   cov = diag(D);  d = sqrt(D)*( inv(L)*V )'*( q - mu );
%
%  Basic S&O algorithm follows
%
   kmx = 2000; lambda = 0; covmx = max( cov ); 
   np = 0; rsqrd = r*r; A = 1;
   for j = 1 : n
     if cov(j) > 1e-10*covmx 
       np = np + 1; gam(np) = 1; 
       alpha(np) = cov(j); A = A/alpha(np); 
       bsqrd(np) = d(j)^2/alpha(np);
       lambda = lambda + bsqrd(np);
     else
       rsqrd = rsqrd - d(j)^2; 
     end
   end
   if ( rsqrd <= 0 ) 
      mvnval = 0; 
   else 
      covmn = min( alpha(1:np) );
      bet = 29*covmn/32; tbeta = rsqrd/bet; A = sqrt(A*bet^np);
      for j = 1 : np
         betalph(j) = bet/alpha(j);
      end
%
%     c(1), c(2), ..., are S&O's c_0, c_1, ...;
%     bsqrd(j) is S&O's b_j^2 ; tbeta is S&O's t/beta, np is S&O's n';
%     g(j) is S&O's g_j; gam(j) is S&O's gamma_j^k.  
%
      csum = A*exp( -lambda/2 ); c(1) = csum; lgb = log( tbeta );  
      if mod( np, 2 )  == 1 
         F = erf( sqrt(tbeta/2) ); lgf = - tbeta/2 - log(2*pi*tbeta)/2;
         for nc = 3:2:np, lgf = lgb + lgf - log(nc-2); F = F - 2*exp(lgf); end
      else  
         F = 1 - exp( -tbeta/2 );  lgf = - tbeta/2 - log(2); 
         for nc = 4:2:np, lgf = lgb + lgf - log(nc-2); F = F - 2*exp(lgf); end
      end
      % equivalent F = gammainc( tbeta/2 , np/2 );  
      mvnval = c(1)*F; 
%
%     for-loop computes S&O series  
%
      for k = 1 : kmx 
         gbsum = 0; 
         for j = 1 : np
            gbsum = gbsum + k*bsqrd(j)*gam(j)*betalph(j);
            gam(j) = gam(j)*( 1 - betalph(j) );
            gbsum = gbsum + gam(j);
         end
         g(k) = gbsum/2; cgsum = 0;
         for j = 0 : k - 1
            cgsum = cgsum + g(k-j)*c(j+1);
         end
         c(k+1) = cgsum/k; csum = csum + c(k+1);
         lgf = lgb + lgf - log(np+2*k-2); F = F - 2*exp(lgf);
         % equivalent F = gammainc( tbeta/2, np/2 + k ); 
         mvnval = mvnval + c(k+1)*F; 
         if ( 1 - csum )*F < re*mvnval, break; end
      end
   end
%
% End function mvnlps
%
