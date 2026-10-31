function [ p, dp, ed ] = qsimvndv( m, r, a, b, t )
%
%  [ p dp e ] = qsimvndv( m, r, a, b, t )
%    uses a randomized quasi-random rule with m points to estimate an
%    MVN probability and t-derivative, for positive definite covariance 
%     matrix r, with lower integration limit column vector t*a 
%     and upper integration limit column vector t*b, and t > 0. 
%   Probability p is output with derivative and error estimates dp and ed.
%    Example usage:
%     >> r = [4 3 2 1;3 5 -1 1;2 -1 4 2;1 1 2 5];
%     >> a = -inf*[1 1 1 1 ]'; b = [ 1 2 3 4 ]';
%     >> [ p d e ] = qsimvndv( 5000, r, a, b, t ); disp([ p d e ])
%
%   This function uses an algorithm given in the paper
%      "Numerical Computation of Multivariate Normal Probabilities", in
%      J. of Computational and Graphical Stat., 1(1992), pp. 141-149, by
%          Alan Genz, WSU Math, PO Box 643113, Pullman, WA 99164-3113
%          Email : alangenz@wsu.edu
%  The primary references for the numerical integration are 
%   "On a Number-Theoretical Integration Method"
%   H. Niederreiter, Aequationes Mathematicae, 8(1972), pp. 304-11, and
%   "Randomization of Number Theoretic Methods for Multiple Integration"
%    R. Cranley and T.N.L. Patterson, SIAM J Numer Anal, 13(1976), pp. 904-14.
%
%   Alan Genz is the author of this function and following Matlab functions.
%
% Initialization
%
[n, n] = size(r); [ ch as bs ] = chlrdr( r, t*a, t*b ); 
ct = ch(1,1); ai = as(1); bi = bs(1);
if ai > -9*ct, if ai < 9*ct, c = phi(ai/ct); else, c=1; end, else c=0; end  
if bi > -9*ct, if bi < 9*ct, d = phi(bi/ct); else, d=1; end, else d=0; end
ci = c; dci = d - ci; p = 0; ed = 0; dp = 0;
ns = 12; nv = fix( max( [ m/ns 1 ] ) ); 
%q = [ 2.^( [1:n-1]/n) 1/nv]'; % Niederreiter point set generators
ps = sqrt(primes(5*n*log(n+1)/4)); q = [ps(1:n-1) 1/nv]'; % Richtmyer generators
%
% Randomization loop for ns samples
%
for i = 1 : ns
   % periodizing transformation 
   xx(:,1:nv) = abs( 2*mod( q*[1:nv] + rand(n,1)*ones(1,nv), 1 ) - 1 );
   [ vp vd ] = mvndnsd( n, nv, ch, ci, dci, xx, as, bs ); 
   d = ( mean(vp) - p )/i;  p = p + d; 
   d = ( mean(vd)/t - dp )/i; dp = dp + d; ed = ( i - 2 )*ed/i + d^2;  
end 
%
ed = 3*sqrt(ed); % error estimate is 3 x standard error with ns samples.
return
%
% end qsimvn
%
function [ p dp ] = mvndnsd( n, nv, ch, ci, dci, x, a, b )
%
%  Transformed integrand for computation of MVN probabilities. 
%
y = zeros(n,nv); c = ci*ones(1,nv); dc = dci*ones(1,nv); p = dc; 
for i = 2 : n
   y(i-1,:) = phinv( c + x(i-1,:).*dc ); s = ch(i,1:i-1)*y(1:i-1,:); 
   ct = ch(i,i)*ones(1,nv); ai = a(i) - s; bi = b(i) - s;
   c = ones( 1, nv ); d = c; 
   c( find( ai <= -9*ct ) ) = 0; d( find( bi <= -9*ct ) ) = 0; 
   tstl = find( ai > -9*ct & ai < 9*ct ); c(tstl) = phi( ai(tstl)./ct(tstl) ); 
   tstl = find( bi > -9*ct & bi < 9*ct ); d(tstl) = phi( bi(tstl)./ct(tstl) ); 
   dc = d - c; p = p.*dc; 
end, y(n,:) = phinv( c + x(n,:).*dc ); dp = ( n - sum(y.*y) ).*p; 

%
% end mvndns
%
function [ c, ap, bp ] = chlrdr( R, a, b )
%
%  Computes permuted lower Cholesky factor c for R which may be singular, 
%   also permuting integration limit vectors a and b.
%
ep = 1e-10; % singularity tolerance;
%
[n,n] = size(R); c = R; ap = a; bp = b; d = sqrt(max(diag(c),0));
for i = 1 :  n
  if d(i) > 0
    c(:,i) = c(:,i)/d(i); c(i,:) = c(i,:)/d(i); 
    ap(i) = ap(i)/d(i); bp(i) = bp(i)/d(i);
  end
end
y = zeros(n,1); sqtp = sqrt(2*pi);
for k = 1 : n
   im = k; ckk = 0; dem = 1; s = 0; 
   for i = k : n 
       if c(i,i) > eps
          cii = sqrt( max( [c(i,i) 0] ) ); 
          if i > 1, s = c(i,1:k-1)*y(1:k-1); end
          ai = ( ap(i)-s )/cii; bi = ( bp(i)-s )/cii; de = phi(bi) - phi(ai);
          if de <= dem, ckk = cii; dem = de; am = ai; bm = bi; im = i; end
       end
   end
   if im > k
      tv = ap(im); ap(im) = ap(k); ap(k) = tv;
      tv = bp(im); bp(im) = bp(k); bp(k) = tv;
      c(im,im) = c(k,k); 
      t = c(im,1:k-1); c(im,1:k-1) = c(k,1:k-1); c(k,1:k-1) = t; 
      t = c(im+1:n,im); c(im+1:n,im) = c(im+1:n,k); c(im+1:n,k) = t; 
      t = c(k+1:im-1,k); c(k+1:im-1,k) = c(im,k+1:im-1)'; c(im,k+1:im-1) = t'; 
   end, c(k,k+1:n) = 0;
   if ckk > ep*k^2
      c(k,k) = ckk; 
      for i = k+1 : n
         c(i,k) = c(i,k)/ckk; c(i,k+1:i) = c(i,k+1:i) - c(i,k)*c(k+1:i,k)';
      end
      if abs(dem) > ep 
	y(k) = ( exp( -am^2/2 ) - exp( -bm^2/2 ) )/( sqtp*dem ); 
      else
	if am < -10
	  y(k) = bm;
	elseif bm > 10
	  y(k) = am;
	else
	  y(k) = ( am + bm )/2;
	end
      end
   else
      c(k:n,k) = 0; y(k) = 0;
   end
end
return
%
% end chlrdr
%
function p = phi(z)
%
%  Standard statistical normal distribution
%
  p = erfc( -z/sqrt(2) )/2;
return
%
% end phi
%
function z = phinv(w)
%
%  Standard statistical inverse normal distribution
%
  z = -sqrt(2)*erfcinv( 2*w );
return
%
% end phinv

