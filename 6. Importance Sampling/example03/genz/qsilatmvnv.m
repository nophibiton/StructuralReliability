function [ p, e ] = qsilatmvnv( m, r, a, b )
%  
%  [ P E ] = QSILATMVNV( M, R, A, B )
%    uses a randomized lattice rule with m points to estimate an
%    MVN probability for positive definite covariance matrix r,
%     with lower integration limit column vector a and upper
%     integration limit column vector b. 
%   Probability p is output with error estimate e.
%    Example use:
%     >> r = [4 3 2 1;3 5 -1 1;2 -1 4 2;1 1 2 5];
%     >> a = -inf*[1 1 1 1 ]'; b = [ 1 2 3 4 ]';
%     >> [ p e ] = qsilatmvnv( 5000, r, a, b ); disp([ p e ])
%
%   This function uses a modification of an algorithm in the paper
%      "Numerical Computation of Multivariate Normal Probabilities", 
%      A. Genz, J. of Comp. Graph. Stat., 1(1992), pp. 141-149.
%
%   The primary references for the numerical integration are 
%    "Randomization of Number Theoretic Methods for Multiple Integration",
%     R. Cranley & T.N.L. Patterson, SIAM J Numer Anal, 13(1976), pp. 904-14.
%   and  
%    "Fast Component-by-Component Construction, a Reprise for Different 
%     Kernels", D. Nuyens & R. Cools. In H. Niederreiter and D. Talay,
%     editors, Monte-Carlo and Quasi-Monte Carlo Methods 2004, 
%     Springer-Verlag, 2006, 371-385.
%
%   Alan Genz is the author of this function and following Matlab functions.
%          Alan Genz, WSU Math, PO Box 643113, Pullman, WA 99164-3113
%          Email : alangenz@wsu.edu
%
% Initialization
%
[n, n] = size(r); [ ch as bs ] = chlrdl( r, a, b ); ct = ch(1,1); 
t = as(1); c = ( 1 + sign(t) )/2; if abs(t) < 9*ct, c = phi(t/ct); end
t = bs(1); d = ( 1 + sign(t) )/2; if abs(t) < 9*ct, d = phi(t/ct); end 
ns = 10; [ q pr ] = fstrnk( m/ns, n-1 ); q = q/pr; y = zeros(n-1,pr); 
ci = c*ones(1,pr); dci = d - ci; p = 0; e = 0; 
%
% Randomization loop for ns samples 
%
for S = 1 : ns, c = ci; dc = dci; vp = dc;
  %
  % Compute randomized MVN integrand with tent periodization
  %
  for i = 2 : n, xi = abs( 2*mod( q(i-1)*[1:pr] + rand, 1 ) - 1 );  
    y(i-1,:) = phinv( c + xi.*dc ); ct = ch(i,i); 
    s = ch(i,1:i-1)*y(1:i-1,:); ai = as(i) - s; bi = bs(i) - s; 
    c = ( 1 + sign(ai) )/2; d = ( 1 + sign(bi) )/2;
    tl = find( abs(ai) < 9*ct ); c(tl) = phi(ai(tl)/ct);
    tl = find( abs(bi) < 9*ct ); d(tl) = phi(bi(tl)/ct);
    dc = d - c; vp = vp.*dc; 
  end, d = ( mean(vp) - p )/S; p = p + d; e = ( S - 2 )*e/S + d^2; 
end, e = 3*sqrt(e); % error estimate is 3 x standard error with ns samples.
%
% end qsilatmvnv
%
function [ c, ap, bp ] = chlrdl( r, a, b )
%
%  Computes scaled permuted lower Cholesky factor c for r which may be 
%  singular, also scaling and permuting integration limit vectors a and b.
%
ep = 1e-10; % singularity tolerance;
%
[n,n] = size(r); c = r; ap = a; bp = b; dc = sqrt(max(diag(c),0)); 
for i = 1 : n, d = dc(i);
  if d > 0, ap(i) = ap(i)/d; bp(i) = bp(i)/d;
    c(:,i) = c(:,i)/d; c(i,:) = c(i,:)/d; 
  end
end, y = zeros(n,1); sqtp = sqrt(2*pi);
for k = 1 : n, epk = ep*k^2; im = k; ck = 0; vm = 1 + epk; 
  for i = k : n 
    if c(i,i) > ep, ci = sqrt(c(i,i)); s = c(i,1:k-1)*y(1:k-1); 
      ai = ( ap(i) - s )/ci; bi = ( bp(i) - s )/ci; 
      dna = 0; dsa = 0; dnb = 0; dsb = 1;
      if ai > -9, dna = exp(-ai^2/2)/sqtp; dsa = phi(ai); end
      if bi <  9, dnb = exp(-bi^2/2)/sqtp; dsb = phi(bi); end, p = dsb - dsa;
      if p > epk, mn = dna - dnb; v = ai*dna - bi*dnb; 
        if      ai < -9, mn = -dnb; v = -bi*dnb;
        elseif  bi >  9, mn =  dna; v =  ai*dna; 
        end, mn = mn/p; v = 1 + v/p - mn^2;
      else, v = 0; mn = ( ai + bi )/2; 
	if ai < -9, mn = bi; elseif bi > 9, mn = ai; end
      end, if v < vm, im = i; vm = v; y(k) = mn; ck = ci; end
    end
  end
  if im > k, c(im,im) = c(k,k); 
    t = c(im,1:k-1); c(im,1:k-1) = c(k,1:k-1); c(k,1:k-1) = t; 
    t = c(im+1:n,im); c(im+1:n,im) = c(im+1:n,k); c(im+1:n,k) = t; 
    t = c(k+1:im-1,k); c(k+1:im-1,k) = c(im,k+1:im-1)'; c(im,k+1:im-1) = t'; 
    tv = ap(im); ap(im) = ap(k); ap(k) = tv;
    tv = bp(im); bp(im) = bp(k); bp(k) = tv;
  end, c(k,k+1:n) = 0; 
  if ck < epk, c(k:n,k) = 0; y(k) = 0;
  else, c(k,k) = ck; c(k+1:n,k) = c(k+1:n,k)/ck; 
    for i = k+1 : n, c(i,k+1:i) = c(i,k+1:i) - c(i,k)*c(k+1:i,k)'; end
  end
end
%
% end chlrdr
%
function p = phi(z)
%
%  Standard statistical normal distribution
%
  p = erfc( -z/sqrt(2) )/2;
%
% end phi
%
function z = phinv(w)
%
%  Standard statistical inverse normal distribution
%
  z = -sqrt(2)*erfcinv( 2*w );
%
% end phinv
%
function [ z, n ] = fstrnk( ni, sm, om, gm, bt )
% 
% Reference: 
%   "Fast Component-by-Component Construction, a Reprise for Different 
%     Kernels", D. Nuyens and R. Cools. In H. Niederreiter and D. Talay,
%     editors, Monte-Carlo and Quasi-Monte Carlo Methods 2004, 
%     Springer-Verlag, 2006, 371-385.
% Modifications to original by A. Genz, 05/07
% Typical Use:  
%  om = inline('x.^2-x+1/6'); n = 99991; s = 100; gam = 0.9.^[1:s];
%  z = fastrank( n, s, om, gam, 1 + gam/3 ); disp([z'; e])
%
n = fix(ni); if ~isprime(n), pt = primes(n); n = pt(length(pt)); end
if nargin < 3, om = inline('x.^2-x+1/6'); 
  bt = ones(1,sm); gm = [ 1 (4/5).^[0:sm-2] ]; 
end, q = 1; w = 1; z = [1:sm]'; m = ( n - 1 )/2; g = prmrot(n); 
perm = [1:m]; for j = 1 : m-1, perm(j+1) = mod( g*perm(j), n ); end
perm = min( n - perm, perm ); c = om(perm/n); fc = fft(c); 
for s = 2 : sm, q = q.*( bt(s-1) + gm(s-1)*c([w:-1:1 m:-1:w+1]) );
  [ es w ] = min( real( ifft( fc.*fft(q) ) ) ); z(s) = perm(w); 
end
%
% end fstrnk
%
function r = prmrot(pin)
%
% find primitive root for prime p, might fail for large primes (p > 32e7)
%
p = pin; if ~isprime(p), pt = primes(p); p = pt(length(pt)); end
pm = p - 1; fp = unique(factor(pm)); n = length(fp); r = 2; k = 1;
while k <= n; d = pm/fp(k); rd = r;
  for i = 2 : d, rd = mod( rd*r, p ); end % computes r^d mod p
  if rd == 1, r = r + 1; k = 0; end, k = k + 1;
end    
%
% prmrot
%





