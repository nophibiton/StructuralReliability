function [ p, e ] = qscmvt( m, nu, sigma, a, cn, b )
%
%  QSCMVT
%    [ p, e ] = qscmvt( m, nu, sigma, a, cn, b )
%    Qsimvt estimates a multivariate t-probability using an
%        m point randomized quasi-random integration method  
%     for a t-distribution with
%        nu degrees-of-freedom, and
%        sigma, a positive semi-definite covariance matrix.       
%     Note: if nu <= 0, a multivariate normal probability is computed.
%    The probability is computed over an integration region specified 
%     by  constraints a < cn*x < b. If sigma is n x n and 
%        cn is k x n, then
%        a and b must be column k-vectors.
%    Output for qscmvt is
%        p, the estimated probability, along with 
%        e, an absolute error estimate for p.
%
%    This function uses an algorithm given in the paper
%     "Comparison of Methods for the Numerical Computation of Multivariate 
%     t Probabilities", J. of Comput. Graph. Stat., 11(2002), pp. 950-971, 
%        by Alan Genz and Frank Bretz. 
%
%   The primary references for the numerical integration are 
%    "On a Number-Theoretical Integration Method"
%    H. Niederreiter, Aequationes Mathematicae, 8(1972), pp. 304-11, and
%    "Randomization of Number Theoretic Methods for Multiple Integration"
%     R. Cranley & T.N.L. Patterson, SIAM J Numer Anal, 13(1976), pp. 904-14.
%
%   Alan Genz is the author of this function and following Matlab functions.
%          Alan Genz, WSU Math, PO Box 643113, Pullman, WA 99164-3113
%          Email : AlanGenz@wsu.edu
%
%
%   Copyright (C) 2005, Alan Genz,  All rights reserved.               
%
%   Redistribution and use in source and binary forms, with or without
%   modification, are permitted provided the following conditions are met:
%     1. Redistributions of source code must retain the above copyright
%        notice, this list of conditions and the following disclaimer.
%     2. Redistributions in binary form must reproduce the above copyright
%        notice, this list of conditions and the following disclaimer in the
%        documentation and/or other materials provided with the distribution.
%     3. The contributor name(s) may not be used to endorse or promote 
%        products derived from this software without specific prior written 
%        permission.
%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
%   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
%   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
%   FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
%   COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
%   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
%   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS 
%   OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
%   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
%   TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
%   USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
[ as ch bs clg n ] = chlsrt( sigma, a, cn, b ); p = 0; e = 0;
ns = 8; nv = fix( max( m/( 2*ns ), 1 ) ); q = 2.^( [1:n]'/(n+1)) ; 
for i = 1 : ns
  vi = 0; xr = rand( n, 1 ); 
  for  j = 1 : nv
    x = abs( 2*mod( j*q + xr, 1 ) - 1 ); 
    if nu > 0, s = sqrt(2*gaminv( x(n) ,nu/2)/nu); else, s = 1; end
    vp =   mvtdnt( n, s*as, ch, s*bs, clg,  x );
    if nu > 0, s = sqrt(2*gaminv(1-x(n),nu/2)/nu); else, s = 1; end
    vp = ( mvtdnt( n, s*as, ch, s*bs, clg, 1-x ) + vp )/2;
    vi = vi + ( vp - vi )/j; 
  end   
  d = ( vi - p )/i; p = p + d; 
  if abs(d) > 0 
    e = abs(d)*sqrt( 1 + ( e/d )^2*( i - 2 )/i );
  else
    if i > 1, e = e*sqrt( ( i - 2 )/i ); end
  end
end
e = 3*e; % error estimate is 3 times standard error with ns samples.
% end qscmvt

function p = mvtdnt( n, a, ch, b, clg, x )
%
%  Transformed integrand for computation of MVT probabilities. 
%
c = phi(a(1)); dc = phi(b(1)) - c; p = dc; 
y = zeros(n-1,1); li = 2; lf = 1;
for i = 2 : n
  y(i-1) = phinv( c + x(i-1)*dc ); lf = lf + clg(i);
  if lf < li, c = 0; dc = 1;
  else, s = ch(li:lf,1:i-1)*y(1:i-1); 
    ai = max( max( a(li:lf) - s ), -9 ); 
    bi = max( ai, min( min( b(li:lf) - s ),  9 ) );
    c = phi(ai); dc = phi(bi) - c; p = p*dc; 
  end, li = li + clg(i);
end 
% end mvtdnt

function [ ap, ch, bp, clg, np ] = chlsrt( r, a, cn, b )
%
%  Computes permuted lower Cholesky factor ch for covariance r which 
%   may be singular, combined with contraints a < cn*x < b, to
%   form revised lower triangular constraint set ap < ch*x < bp; 
%   clg contains information about structure of ch: clg(1) rows for 
%   ch with 1 nonzero, ..., clg(np) rows with np nonzeros.
%
ep = 1e-10; % singularity tolerance;
%
[n n] = size(r); [m n] = size(cn); ch = cn; np = 0;
ap = a; bp = b; y = zeros(n,1); sqtp = sqrt(2*pi);
c = r; d = sqrt(max(diag(c),0));
for i = 1 : n, di = d(i);
  if di > 0
    c(:,i) = c(:,i)/di; c(i,:) = c(i,:)/di; ch(:,i) = ch(:,i)*di;
  end
end
%
% determine (with pivoting) Cholesky factor for r 
%  and form revised constraint matrix ch
%
for i = 1 : n
  clg(i) = 0; epi = ep*i^2; j = i; 
  for l = i+1 : n, if c(l,l) > c(j,j), j = l; end, end
  if j > i
    t = c(i,i); c(i,i) = c(j,j); c(j,j) = t;
    t = c(i,1:i-1); c(i,1:i-1) = c(j,1:i-1); c(j,1:i-1) = t;
    t = c(i+1:j-1,i); c(i+1:j-1,i) = c(j,i+1:j-1)'; c(j,i+1:j-1) = t';
    t = c(j+1:n,i); c(j+1:n,i) = c(j+1:n,j); c(j+1:n,j) = t;
    t = ch(:,i); ch(:,i) = ch(:,j); ch(:,j) = t;
  end
  if c(i,i) < epi, break, end, cvd = sqrt( c(i,i) ); c(i,i) = cvd;
  for l = i+1 : n
    c(l,i) = c(l,i)/cvd; c(l,i+1:l) = c(l,i+1:l) - c(l,i)*c(i+1:l,i)';
  end
  ch(:,i) = ch(:,i:n)*c(i:n,i); np = np + 1;
end
%
% use right reflectors to reduce ch to lower triangular
%
for i = 1 : min( np-1, m )
  epi = ep*i*i; vm = 1; lm = i;
  %
  % permute rows so that smallest variance variables are first.
  %
  for l = i : m    
    v = ch(l,1:np); s = v(1:i-1)*y(1:i-1); 
    ss = max( sqrt( sum( v(i:np).^2 ) ), epi ); 
    al = ( ap(l) - s )/ss; bl = ( bp(l) - s )/ss; 
    dna = 0; dsa = 0; dnb = 0; dsb = 1;
    if al > -9, dna = exp(-al*al/2)/sqtp; dsa = phi(al); else, al = bl; end
    if bl <  9, dnb = exp(-bl*bl/2)/sqtp; dsb = phi(bl); else, bl = al; end
    if dsb - dsa > epi, ds = dsb - dsa; 
      mn = ( dna - dnb )/ds; vr = 1 + ( al*dna - bl*dnb )/ds - mn^2; 
    else, mn = ( al + bl )/2; vr = 0;
    end, if vr <= vm, lm = l; vm = vr; y(i) = mn; end
  end
  v = ch(lm,1:np); 
  if lm > i 
    ch(lm,1:np) = ch(i,1:np); ch(i,1:np) = v;
    tl = ap(i); ap(i) = ap(lm); ap(lm) = tl;
    tl = bp(i); bp(i) = bp(lm); bp(lm) = tl;
  end
  ch(i,i+1:np) = 0; ss = sum( v(i+1:np).^2 );
  if ( ss > epi )
    ss = sqrt( ss + v(i)^2 ); if v(i) < 0, ss = -ss; end
    ch(i,i) = -ss; v(i) = v(i) + ss; vt = v(i:np)'/( ss*v(i) );
    ch(i+1:m,i:np) = ch(i+1:m,i:np) - ch(i+1:m,i:np)*vt*v(i:np);
  end
end
%
% scale and sort constraints
%
for i = 1 : m
  v(1:np) = ch(i,1:np); clm(i) = min(i,np); 
  jm = 1; for j = 1 : clm(i), if abs(v(j)) > ep*j*j, jm = j; end, end 
  if jm < np, v(jm+1:np) = 0; end, clg(jm) = clg(jm) + 1;
  at = ap(i); bt = bp(i); j = i;
  for l = i-1 : -1 : 1
    if jm >= clm(l), break, end
    ch(l+1,1:np) = ch(l,1:np); j = l;
    ap(l+1) = ap(l); bp(l+1) = bp(l); clm(l+1) = clm(l);
  end
  clm(j) = jm; vjm = v(jm); ch(j,1:np) = v/vjm;
  ap(j) = at/vjm; bp(j) = bt/vjm;
  if vjm < 0, tl = ap(j); ap(j) = bp(j); bp(j) = tl; end
end
j = 0; for i = 1 : np, if clg(i) > 0, j = i; end, end, np = j;
%
% combine constraints for first variable
%
if clg(1) > 1 
  ap(1) = max( ap(1:clg(1)) ); bp(1) = max( ap(1), min( bp(1:clg(1)) ) ); 
  ap(2:m-clg(1)+1) = ap(clg(1)+1:m); bp(2:m-clg(1)+1) = bp(clg(1)+1:m);
  ch(2:m-clg(1)+1,:) = ch(clg(1)+1:m,:); 
  clg(1) = 1;
end
%
% end chlsrt
%

function p = phi(z)
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
