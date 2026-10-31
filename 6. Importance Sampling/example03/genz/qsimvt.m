function [ p, e ] = qsimvt( m, nu, sigma, a, b )
%
%  QSIMVT
%    Qsimvt estimates a multivariate t-probability using an
%        m point randomized quasi-random integration method  
%    for a t-distribution with
%        nu degrees-of-freedom, and
%        sigma, a positive definite covariance matrix.
%    The probability is computed over a hyper-rectangle specified with
%        a, lower integration limits and 
%        b, upper integration limits. 
%      If sigma is nxn then a and b must be n-vectors.
%    Output for qsimvt is
%        p, the estimated probability, along with 
%        e, an absolute error estimate for p.
%    Example:
%     >> r = [4 3 2 1;3 5 -1 1;2 -1 4 2;1 1 2 5];
%     >> a = -inf*[1 1 1 1 ]'; b = [ 1 2 3 4 ]';
%     >> [ p e ] = qsimvt( 5000, 5, r, a, b ); disp([ p e ])

%
%   This function uses an algorithm given in the paper
%      "Comparison of Methods for the Numerical Computation of 
%       Multivariate t Probabilities", in
%      J. of Computational and Graphical Stat., 11(2002), pp. 950-971, by
%          Alan Genz and Frank Bretz 
%
%   The primary references for the numerical integration are 
%    "On a Number-Theoretical Integration Method"
%    H. Niederreiter, Aequationes Mathematicae, 8(1972), pp. 304-11.
%    and
%    "Randomization of Number Theoretic Methods for Multiple Integration"
%     R. Cranley & T.N.L. Patterson, SIAM J Numer Anal, 13(1976), pp. 904-14.
%
%   Alan Genz is the author of qsimvt and following Matlab functions.
%          Alan Genz
%          Department of Mathematics
%          Washington State University 
%          Pullman, WA 99164-3113
%          Email : AlanGenz@wsu.edu
%

% Initialization
%
[n, n] = size(sigma); 
sn = max([1 sqrt(nu)]); [ ch as bs ] = chlrps( sigma, a/sn, b/sn ); 
ns = 8; nv = max( [ m/( 2*ns ) 1 ] ); p = 0; e = 0;
% q = 2.^( [1:n]'/(n+1)) ; % Niederreiter point set generators
ps = sqrt(primes(5*(n+1)*log(n+2)/4)); q = ps(1:n)'; % Richtmyer generators
%
% Randomization loop for ns samples
%
for i = 1 : ns
  vi = 0; xr = rand( n, 1 ); 
  %
  % Loop for 2*nv quasirandom points
  %  
  for  j = 1 : nv
    x = abs( 2*mod( j*q + xr, 1 ) - 1 );  % periodizing transformation
    if nu > 0, s = sqrt( 2*gaminv(  x(n),  nu/2 ) ); else s = 1; end 
    at = s*as; bt = s*bs; 
    vp =   mvtdns( n, ch, at, bt,  x );
    if nu > 0, s = sqrt( 2*gaminv( 1-x(n), nu/2 ) ); else s = 1; end 
    at = s*as; bt = s*bs; 
    vp = ( mvtdns( n, ch, at, bt, 1-x ) + vp )/2;
    vi = vi + ( vp - vi )/j; 
  end   
  %
  d = ( vi - p )/i; p = p + d; 
  if abs(d) > 0 
    e = abs(d)*sqrt( 1 + ( e/d )^2*( i - 2 )/i );
  else
    if i > 1, e = e*sqrt( ( i - 2 )/i ); end
  end  
end
%
e = 3*e; % error estimate is 3 times standard error with ns samples.
% end qsimvt
%
function p = mvtdns( n, ch, a, b, x )
%
%  Transformed integrand for computation of MVT probabilities. 
%
y = zeros(n-1,1); cn = 37.5; s = 0; p = 1; 
for i = 1 : n
  if i > 1, y(i-1) = phinv( c + x(i-1)*dc ); s = ch(i,1:i-1)*y(1:i-1); end
  ai = a(i) - s; bi = b(i) - s;
  if abs(ai) < cn, c = phid(ai); else, c = ( 1 + sign(ai) )/2; end
  if abs(bi) < cn, d = phid(bi); else, d = ( 1 + sign(bi) )/2; end
  dc = d - c; p = p*dc; 
end 
%
% end mvtdns
%
function [ c, ap, bp ] = chlrps( R, a, b )
%
%  Computes permuted and scaled lower Cholesky factor c for R which may be 
%   singular, also permuting and scaling integration limit vectors a and b.
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
          cii = sqrt( max( c(i,i), 0 ) ); 
          if i > 1, s = c(i,1:k-1)*y(1:k-1); end
          ai = (ap(i)-s)/cii; bi = (bp(i)-s)/cii; de = phid(bi) - phid(ai);
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
   end
   if ckk > ep*k^2
      c(k,k) = ckk; c(k,k+1:n) = 0;
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
      c(k,1:k) = c(k,1:k)/ckk; ap(k) = ap(k)/ckk; bp(k) = bp(k)/ckk; 
   else
      c(k:n,k) = 0; y(k) = ( ap(k) + bp(k) )/2;
   end
end
%
function p = phid(z)
%
%  Standard statistical normal distribution
%
  p = erfc( -z/sqrt(2) )/2; 
%
function z = phinv(w)
%
%  Standard statistical inverse normal distribution
%
  z = -sqrt(2)*erfcinv( 2*w );
%
