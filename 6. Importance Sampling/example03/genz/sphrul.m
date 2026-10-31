function [ intval, intcls ] = sphrul( s, d, f )
%
%  SPHRUL computes a fully symmetric rule approximation to an 
%   integral over the surface of the unit sphere: norm(x) = 1. 
%**************  parameters for sphrul  ********************************
%***** input parameters
%  s       integer number of variables
%  d       integer degree parameter
%  f       user defined real function integrand f(x).
%****** output parameters
%  intval  intval will be an approximation of polynomial degree 2d + 1.
%  intcls  integer total number of f values needed for intval
%****** example usage
% >> ft = inline('x(1)^2*x(2)^4*cos(x(4))*exp(-x(3))','x');
% >> for d = 2 : 6, [ int pts ] = sphrul(4,d,ft); disp([d int pts]), end
% 
%     Sphrul uses a transformation of Sylvester's interpolation formula
%     for simplices, as described in "Fully Symmetric Interpolatory
%     Rules for Multiple Integrals over Hyper-Spherical Surfaces",
%       J. Comp. Appl. Math. 157 (2003), pp. 187-195
%   Software and paper written by 
%     Alan Genz
%     Department of Mathematics
%     Washington State University
%     Pullman, Washington 99164-3113  USA
%     email: alangenz@wsu.edu
%
%***********************************************************************
%
%***  begin loop for each d
%      for each d find all distinct partitions m with mod(m) = d
%
  m(1) = d; intval = 0; intcls = 0; for i = 2 : s, m(i) = 0; end
  while m(1) <= d 
    %     
    %***  calculate weight for partition m and fully symmetric sums
    %
    wt = weight( s, m, d ); 
    if abs(wt) > 1d-10, [ fs sumcls ] = fulsum( s, m, d, f ); 
      intval = intval + wt*fs; intcls = intcls + sumcls;
    end
    m = nxprtr( s, m );
  end
  intval = srfsph(s)*intval;
%
%  end function sphrul
%
function v = srfsph( s )
%
% compute surface content for s-sphere
%
   v = 2; if mod(s,2) == 0, v = v*pi; end
   for i = s-2 : -2 : 1, v = 2*v*pi/i; end
%
%*
function ws = weight( s, m, d )
%
%***  calculate weight for partition m
%
  ws = 0; for i = 1: d, k(i) = 0; end  
  while 1
    is = 0; ks = 0; mc = 1; kn = 1; pd = 1;
    for i = 1 : d, is = is + 1;
      if k(i) == 0, pd = pd*( 1 - is );
      else, pd = d*kn*pd/( s + ks ); ks = ks + 2; kn = kn + 2;
      end 
      if is == m(mc), mc = mc + 1; is = 0; kn = 1; end
    end
    ws = ws + pd;
    for i = 1 : d
      k(i) = k(i) + 1; if k(i) < 2, break, end, k(i) = 0;
    end
    if i == d & k(d) == 0, break, end
  end
  for i = 1 : s, for j = 1 : m(i), ws = ws/( m(i) - j + 1 ); end, end
%
function [ fs, sumcls ] = fulsum( s, m, d, f )
%
%***  compute fully symmetric basic rule sum
%
  sumcls = 0; fs = 0; fin = 0;
  %
  %*******  compute centrally symmetric sum for permutation m
  %
  while fin == 0, x = -sqrt( m/d );
    [ sym, symcls ] = symsum( s, x, f );
    fs = fs + sym; sumcls = sumcls + symcls;
    %     
    %*******  find next distinct permutation of m and loop back
    %          to compute next centrally symmetric sum
    %
    [ m, fin ] = nxperm( s, m );
  end
%     
%*****end loop for permutations of m and associated sums
%
%  end function fulsum
%
function [ sym, sumcls ] = symsum( s, zi, f )
%
%***  function to compute symmetric basic rule sum
%
  sumcls = 0; sym = 0; z = zi;
  %
  %*******  compute centrally symmetric sum for zi
  %
  while 1
    sumcls = sumcls + 1; sym = sym + ( feval( f, z ) - sym )/sumcls;
    for i = 1 : s, z(i) = -z(i); if z(i) > 0, break, end, end
    if z == zi, break, end
  end
%
%  end function symsum
%
%
function [ m, fin ] = nxperm( s, mi )
%
%***  compute the next permutation of mi
%
  m = mi; fin = 0;
  for i = 2 : s
    if m(i-1) > m(i), mi = m(i); ixchng = i - 1;
      if i > 2
	for l = 1 : ixchng/2
	  ml = m(l); imnusl = i - l; m(l) = m(imnusl); m(imnusl) = ml;
          if ml  <= mi, ixchng = ixchng - 1; end
	  if m(l) > mi, lxchng = l; end
	end
	if m(ixchng) <= mi, ixchng = lxchng; end
      end 
      m(i) = m(ixchng); m(ixchng) = mi; return
    end 
  end
  %
  %*** restore original order to m.
  %
  for i = 1 : s/2, mi = m(i); m(i) = m(s-i+1); m(s-i+1) = mi; end
  fin = 1;
%
%  end function nxperm
%
function p = nxprtr( s, pi )
%
%***  compute the next s partition of pi
%
  p = pi; psum = p(1);
  for i = 2 : s, psum = psum + p(i);
    if p(1) <= p(i) + 1, p(i) = 0;
    else
      p(1) = psum - (i-1)*( p(i) + 1 );
      for l = 2 : i, p(l) = p(i) + 1; end
      return
    end
  end
  p(1) = psum + 1;
%
%  end function nxtptr
%