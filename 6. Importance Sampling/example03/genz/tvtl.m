%
%  This file contains functions tvtl (trivariate normal and t), 
%  bvtl (bivariate t), bvnl (bivariate normal), studnt (univariate t), 
%  phid (univariate normal), plus support functions.
%
%  The software is based on work described in the paper
%   "Numerical Computation of Rectangular Bivariate and Trivariate 
%    Normal and t Probabilities", by
%          Alan Genz
%          Department of Mathematics
%          Washington State University
%          Pullman, WA 99164-3113
%          Email : alangenz@wsu.edu
%
function tvt = tvtl( nu, h, r, epsi )
%    
%     A function for computing trivariate normal and t-probabilities.
%     tvtl calculates the probability that x(i) < h(i), for i = 1, 2, 3.     
%   parameters
%    nu    integer degrees of freedom; use nu = 0 for normal cases.
%    h     real array of three upper limits for probability distribution 
%    r     real array of three correlation coefficients, r should 
%          contain the lower left portion of the correlation matrix R. 
%          r should contain the values r21, r31, r32 in that order.
%    epsi  real required absolute accuracy; maximum accuracy for most
%          computations is approximately 1e-14.
% 
%     This function uses algorithms developed from the ideas 
%     described in the papers:
%       R.L. Plackett, Biometrika 41(1954), pp. 351-360.
%       Z. Drezner, Math. Comp. 62(1994), pp. 289-294.
%     with adaptive integration from (0,0,1) to (0,0,r32) to R. 
%
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
  global TVTH1 TVTH2 TVTH3 TVTR23 TVTRUA TVTRUB TVTAR TVTRUC TVTNUC
  epst = max( [ 1e-14 epsi ] );
  pt = pi/2;
  h1 = h(1);
  h2 = h(2);
  h3 = h(3);
  r12 = r(1);
  r13 = r(2);
  r23 = r(3);
  %
  %     Sort R's and check for special cases
  %
  if ( abs(r12) > abs(r13) ) 
    h2 = h3;
    h3 = h(2);
    r12 = r13;
    r13 = r(1);
  end
  if ( abs(r13) > abs(r23) ) 
    h1 = h2;
    h2 = h(1);
    r23 = r13;
    r13 = r(3);
  end
  tvt = 0;
  if ( abs(h1) + abs(h2) + abs(h3) < epst )
    tvt = ( 1 + ( asin(r12) + asin(r13) + asin(r23) )/pt )/8;
  elseif ( nu < 1 & abs(r12) + abs(r13) < epst )
    tvt = phid(h1)*bvtl( nu, h2, h3, r23 );
  elseif ( nu < 1 & abs(r13) + abs(r23) < epst ) 
    tvt = phid(h3)*bvtl( nu, h1, h2, r12 );
  elseif ( nu < 1 & abs(r12) + abs(r23) < epst ) 
    tvt = phid(h2)*bvtl( nu, h1, h3, r13 );
  elseif ( 1 - r23 < epst ) 
    tvt = bvtl( nu, h1, min( [h2 h3] ), r12 );
  elseif ( r23 + 1 < epst ) 
    if  ( h2 > -h3 ) 
      tvt = bvtl( nu, h1, h2, r12 ) - bvtl( nu, h1, -h3, r12 );
    end
  else
    %
    %        Compute singular TVT value
    %
    if ( nu < 1 ) 
      tvt = bvtl( nu, h2, h3, r23 )*phid(h1);
    elseif ( r23 >= 0 ) 
      tvt = bvtl( nu, h1, min( [h2 h3] ), 0 );
    elseif ( h2 > -h3 ) 
      tvt = bvtl( nu, h1, h2, 0 ) - bvtl( nu, h1, -h3, 0 );
    end
    %
    %        Use numerical integration to compute probability;
    %
    %
    TVTH1 = h1; 
    TVTH2 = h2; 
    TVTH3 = h3; 
    TVTR23 = r23;
    TVTRUA = asin( r12 );
    TVTRUB = asin( r13 );
    ar = asin( r23);
    TVTAR = ar;
    TVTRUC = pi*sign( ar )/2 - ar;
    TVTNUC = nu;
    tvt = tvt + adonet( 'tvtmfn', 0, 1, epst )/( 2*pi );
  end
  tvt = max( 0, min( tvt, 1 ) ) ;
  return
%
% end tvtl
%
function f = tvtmfn( x )
%
%     Computes Plackett formula integrands
%
  global TVTH1 TVTH2 TVTH3 TVTR23 TVTRUA TVTRUB TVTAR TVTRUC TVTNUC;
  f = 0;
  [ r12 rr2 ] = sincs( TVTRUA*x );
  [ r13 rr3 ] = sincs( TVTRUB*x );
  if ( abs(TVTRUA) > 0 ) 
    f = f + TVTRUA*pntgnd( TVTNUC,TVTH1,TVTH2,TVTH3,r13,TVTR23,r12,rr2 );
  end
  if ( abs(TVTRUB) > 0 )
    f = f + TVTRUB*pntgnd( TVTNUC,TVTH1,TVTH3,TVTH2,r12,TVTR23,r13,rr3 );
  end
  if ( TVTNUC > 0 ) 
    [ r rr ] = sincs( TVTAR + TVTRUC*x );
    f = f - TVTRUC*pntgnd( TVTNUC, TVTH2, TVTH3, TVTH1, 0, 0, r, rr );
  end    
  return
%
% end tvtmfn
%
   function [sx, cs] = sincs( x )
%
%     Computes sin(x), cos(x)^2, with series approx. for |x| near pi/2
%
      ee = ( pi/2 - abs(x) )^2;
      if ( ee < 5e-5 ) 
         sx = ( 1 - ee*( 1 - ee/12 )/2 )*sign( x );
         cs = ee*( 1 - ee*( 1 - 2*ee/15 )/3 );
      else
         sx = sin(x);
         cs = 1 - sx*sx;
      end
   return
%
% end sincs
%
function f = pntgnd( nu, ba, bb, bc, ra, rb, r, rr )
%
%     Computes Plackett formula integrand;
%
  f = 0;
  dt = rr*( rr - ( ra - rb )^2 - 2*ra*rb*( 1 - r ) );
  if ( dt > 0 ) 
    bt = ( bc*rr + ba*( r*rb - ra ) + bb*( r*ra -rb ) )/sqrt(dt) ;
    ft = ( ba - r*bb )^2/rr + bb*bb;
    if ( nu < 1 ) 
      if ( bt > -10 & ft < 100 ) 
	f = exp( -ft/2 );
	if ( bt < 10 ),  f = f*phid(bt); end
      end
    else
      ft = sqrt( 1 + ft/nu );
      f = studnt( nu, bt/ft )/ft^nu;
    end
  end
  return
%
% end pntgnd
%
function st = studnt( nu, t )
%
%     Student t Distribution Function;
%
%                       t
%         studnt = c   i  ( 1 + y*y/nu )^( -(nu+1)/2 ) dy
%                   nu -inf
%  
  if ( nu < 1 ) 
    st = phid( t );
  elseif ( nu == 1 ) 
    st = ( 1 + 2*atan(t)/pi )/2;
  elseif ( nu == 2 ) 
    st = ( 1 + t/sqrt( 2 + t*t ))/2;
  else 
    tt = t*t;
    cssthe = 1/( 1 + tt/nu );
    polyn = 1;
    for j = nu-2 : -2 : 2
      polyn = 1 + ( j - 1 )*cssthe*polyn/j;
    end 
    if ( mod( nu, 2 ) == 1 ) 
      rn = nu;
      ts = t/sqrt(rn);
      st = ( 1 + 2*( atan(ts) + ts*cssthe*polyn )/pi )/2;
    else
      snthe = t/sqrt( nu + tt );
      st = ( 1 + snthe*polyn )/2;
    end
    st = max( [ 0 min( [st 1] )] );
  end
  return
%
% end studnt
%
function p = bvtl( nu, dh, dk, r )
%
%    A function for computing bivariate t probabilities.
%    bvtl calculates the probability that x < dh and y < dk;
%   parameters
%     nu number of degrees of freedom
%     dh 1st upper integration limit
%     dk 2nd upper integration limit
%     r   correlation coefficient
%
%        This function is based on the method described by 
%          Dunnett, C.W. and M. Sobel, (1954),
%          A bivariate generalization of Student's t-distribution
%          with tables for certain special cases,
%          Biometrika 41, pp. 153-169,
%
%       Alan Genz
%       Department of Mathematics
%       Washington State University
%       Pullman, Wa 99164-3113
%       Email : alangenz@wsu.edu
%
  if nu < 1 
    p = bvnl( dh, dk, r );
  elseif ( 1 - r <= eps ) 
    p = studnt( nu, min([ dh dk ]) );
  elseif ( r + 1  <= eps ) 
    if ( dh > -dk )  
      p = studnt( nu, dh ) - studnt( nu, -dk );
    else
      p = 0;
    end 
  else
    tpi = 2*pi; ors = 1 - r*r; hrk = dh - r*dk; krh = dk - r*dh; 
    if abs(hrk) + ors > 0 
      xnhk = hrk^2/( hrk^2 + ors*( nu + dk^2 ) );
      xnkh = krh^2/( krh^2 + ors*( nu + dh^2 ) );
    else
      xnhk = 0; xnkh = 0; 
    end
    hs = sign( dh - r*dk ); ks = sign( dk - r*dh );
    if mod( nu, 2 ) == 0
      bvt = atan2( sqrt(ors), -r )/tpi;
      gmph = dh/sqrt( 16*( nu + dh^2 ) ); gmpk = dk/sqrt( 16*( nu + dk^2 ) ); 
      btnckh = 2*atan2( sqrt( xnkh ), sqrt( 1 - xnkh ) )/pi; 
      btpdkh = 2*sqrt( xnkh*( 1 - xnkh ) )/pi;
      btnchk = 2*atan2( sqrt( xnhk ), sqrt( 1 - xnhk ) )/pi; 
      btpdhk = 2*sqrt( xnhk*( 1 - xnhk ) )/pi;
      for j = 1 : nu/2
	bvt = bvt + gmph*( 1 + ks*btnckh ); 
	bvt = bvt + gmpk*( 1 + hs*btnchk );
	btnckh = btnckh + btpdkh; btpdkh = 2*j*btpdkh*( 1 - xnkh )/(2*j+1); 
	btnchk = btnchk + btpdhk; btpdhk = 2*j*btpdhk*( 1 - xnhk )/(2*j+1); 
	gmph = gmph*( j - 1/2 )/( j*( 1 + dh^2/nu ) );
	gmpk = gmpk*( j - 1/2 )/( j*( 1 + dk^2/nu ) );
      end
    else
      qhrk = sqrt( dh^2 + dk^2 - 2*r*dh*dk + nu*ors ); 
      hkrn = dh*dk + r*nu; hkn = dh*dk - nu; hpk = dh + dk;
      bvt = atan2( -sqrt(nu)*(hkn*qhrk+hpk*hkrn), hkn*hkrn-nu*hpk*qhrk )/tpi; 
      if bvt < -10*eps, bvt = bvt + 1; end
      gmph = dh/( tpi*sqrt(nu)*( 1 + dh^2/nu ) ); 
      gmpk = dk/( tpi*sqrt(nu)*( 1 + dk^2/nu ) ); 
      btnckh = sqrt( xnkh ); btpdkh = btnckh;
      btnchk = sqrt( xnhk ); btpdhk = btnchk; 
      for j = 1 : ( nu - 1 )/2
	bvt = bvt + gmph*( 1 + ks*btnckh ); 
	bvt = bvt + gmpk*( 1 + hs*btnchk );
	btpdkh = (2*j-1)*btpdkh*( 1 - xnkh )/(2*j); btnckh = btnckh + btpdkh; 
	btpdhk = (2*j-1)*btpdhk*( 1 - xnhk )/(2*j); btnchk = btnchk + btpdhk; 
	gmph = gmph*j/( ( j + 1/2 )*( 1 + dh^2/nu ) );
	gmpk = gmpk*j/( ( j + 1/2 )*( 1 + dk^2/nu ) );
      end
    end
    p = bvt;
  end
  return
%
%  end bvtl
%
function p = phid(z)
  p = erfc( -z/sqrt(2) )/2;
  return
%
% end phid
%
function p = bvnl( dh, dk, r )
%
%   a function for computing bivariate normal probabilities.
%   bvnl calculates the probability that x < dh and y < dk. 
%   parameters
%     dh 1st upper integration limit
%     dk 2nd upper integration limit
%     r   correlation coefficient
%
%    this function is based on the method described by 
%        Drezner, Z and G.O. Wesolowsky, (1989),
%        On the computation of the bivariate normal inegral,
%        Journal of Statist. Comput. Simul. 35, pp. 101-107,
%    with major modifications for double precision, for |r| close to 1,
%    and for matlab by Alan Genz - last modifications 7/98.
%       Alan Genz
%       Department of Mathematics
%       Washington State University
%       Pullman, Wa 99164-3113
%       Email : alangenz@wsu.edu
%
  p = bvnu( -dh, -dk, r );
  return
%
%   end bvnl
%
function p = bvnu( dh, dk, r )
%
%  A function for computing bivariate normal probabilities.
%  bvnu calculates the probability that x > dh and y > dk. 
%    parameters  
%      dh 1st lower integration limit
%      dk 2nd lower integration limit
%      r   correlation coefficient
%
%   Author
%       Alan Genz
%       Department of Mathematics
%       Washington State University
%       Pullman, Wa 99164-3113
%       Email : alangenz@wsu.edu
%
%    This function is based on the method described by 
%        Drezner, Z and G.O. Wesolowsky, (1989),
%        On the computation of the bivariate normal inegral,
%        Journal of Statist. Comput. Simul. 35, pp. 101-107,
%    with major modifications for double precision, for |r| close to 1,
%    and for matlab by Alan Genz - last modifications 7/98.
%        Note: to compute the probability that x < dh and y < dk, use 
%              bvnu( -dh, -dk, r ). 
%
  twopi = 2*pi; 
  if ( abs(r) < 0.3 ) 
    ng = 1; lg = 3;
    %       Gauss Legendre points and weights, n =  6
    w(1:3,1) = [0.1713244923791705 0.3607615730481384 0.4679139345726904]';
    x(1:3,1) = [0.9324695142031522 0.6612093864662647 0.2386191860831970]';
  elseif ( abs(r) < 0.75 ) 
    ng = 2; lg = 6;
    %       Gauss Legendre points and weights, n = 12
    w(1:3,2) = [.04717533638651177 0.1069393259953183 0.1600783285433464]';
    w(4:6,2) = [0.2031674267230659 0.2334925365383547 0.2491470458134029]';
    x(1:3,2) = [0.9815606342467191 0.9041172563704750 0.7699026741943050]';
    x(4:6,2) = [0.5873179542866171 0.3678314989981802 0.1252334085114692]';
  else 
    ng = 3; lg = 10;
    %       Gauss Legendre points and weights, n = 20
    w(1:3,3) = [.01761400713915212 .04060142980038694 .06267204833410906]';
    w(4:6,3) = [.08327674157670475 0.1019301198172404 0.1181945319615184]';
    w(7:9,3) = [0.1316886384491766 0.1420961093183821 0.1491729864726037]';
    w(10,3) = 0.1527533871307259;
    x(1:3,3) = [0.9931285991850949 0.9639719272779138 0.9122344282513259]';
    x(4:6,3) = [0.8391169718222188 0.7463319064601508 0.6360536807265150]';
    x(7:9,3) = [0.5108670019508271 0.3737060887154196 0.2277858511416451]';
    x(10,3) = 0.07652652113349733;
  end
  h = dh; k = dk; hk = h*k; bvn = 0;
  if ( abs(r) < 0.925 ) 
    hs = ( h*h + k*k )/2; asr = asin(r);
    for i = 1 : lg
      sn = sin( asr*( 1 - x(i,ng) )/2 );
      bvn = bvn + w(i,ng)*exp( ( sn*hk - hs )/( 1 - sn*sn ) );
      sn = sin( asr*( 1 + x(i,ng) )/2 );
      bvn = bvn + w(i,ng)*exp( ( sn*hk - hs )/( 1 - sn*sn ) );
    end
    bvn = bvn*asr/( 4*pi ) + phid(-h)*phid(-k) ;
  else
    if ( r < 0 ), k = -k; hk = -hk; end
    if ( abs(r) < 1 ) 
      as = ( 1 - r )*( 1 + r ); a = sqrt(as); bs = ( h - k )^2;
      c = ( 4 - hk )/8 ; d = ( 12 - hk )/16; asr = -( bs/as + hk )/2;
      if ( asr > -100 ) 
	bvn = a*exp(asr)*( 1 - c*(bs-as)*(1-d*bs/5)/3 + c*d*as*as/5 );
      end
      if ( -hk < 100 ) 
	b = sqrt(bs); sp = sqrt(twopi)*phid(-b/a);
	bvn = bvn - exp(-hk/2)*sp*b*( 1 - c*bs*(1 - d*bs/5 )/3 );
      end
      a = a/2;
      for i = 1 : lg
	for is = -1 : 2 : 1
	  xs = ( a*(  is*x(i,ng) + 1 ) )^2; rs = sqrt( 1 - xs );
	  asr = -( bs/xs + hk )/2;
	  if ( asr > -100 )
	    sp = ( 1 + c*xs*( 1 + d*xs ) );
	    ep = exp( -hk*( 1 - rs )/( 2*( 1 + rs ) ) )/rs;
	    bvn = bvn + a*w(i,ng)*exp(asr)*( ep - sp );
	  end
	end
      end
      bvn = -bvn/twopi;
    end
    if r > 0 , bvn =  bvn + phid( -max( h, k ) ); end
    if r < 0 , bvn = -bvn + max( 0, phid(-h)-phid(-k) ); end
  end
  p = max( 0, min( 1, bvn ) );
  return
%
%   end bvnu
%
function fin = adonet( f, a, b, tol )
%
%     1-dimensional adaptive integration
%
  nl = 100;
  ai(1) = a;
  bi(1) = b;
  ip = 1;
  im = 1;
  err = 1;
  while ( 4*err > tol & im < nl )
    im = im + 1;
    bi(im) = bi(ip);
    ai(im) = ( ai(ip) + bi(ip) )/2;
    bi(ip) = ai(im);
    [ fi(ip) ei(ip) ] = krnrdt( ai(ip), bi(ip), f );
    [ fi(im) ei(im) ] = krnrdt( ai(im), bi(im), f );
    fin = sum( fi(1:im) );
    err = sqrt( sum( ei(1:im).^2 ) );
    for i = 1 : im, if ei(i) >  ei(ip), ip = i; end, end
  end  
  return
%
function [ resk, err ] = krnrdt( a, b, f )
%
%     kronrod rule
%
  wg0    = 0.2729250867779007e+00;
  wg(1:3) = [ .05566856711617449 0.1255803694649048 0.1862902109277352 ];
  wg(4:5) = [ 0.2331937645919914 0.2628045445102478 ];
%
  xgk(1:3) = [ 0.9963696138895427 0.9782286581460570 0.9416771085780681 ];
  xgk(4:6) = [ 0.8870625997680953 0.8160574566562211 0.7301520055740492 ];
  xgk(7:9) = [ 0.6305995201619651 0.5190961292068118 0.3979441409523776 ];
  xgk(10:11) = [ 0.2695431559523450 0.1361130007993617 ];
%
  wgk0 = 0.1365777947111183e+00;
  wgk(1:3) = [ .00976544104596129 .02715655468210443 .04582937856442671 ];
  wgk(4:6) = [ .06309742475037484 .07866457193222764 .09295309859690074 ];
  wgk(7:9) = [ 0.1058720744813894 0.1167395024610472 0.1251587991003195 ];
  wgk(10:11) = [ 0.1312806842298057 0.1351935727998845 ];
  wid = ( b - a )/2;
  cen = ( b + a )/2;
  fc = feval( f, cen );
  resg = fc*wg0;
  resk = fc*wgk0;
  for j = 1 : 5
    t = wid*xgk(2*j-1); 
    fc = feval( f, cen - t ) + feval( f, cen + t );
    resk = resk + wgk(2*j-1)*fc;
    t = wid*xgk(2*j); 
    fc = feval( f, cen - t ) + feval( f, cen + t );
    resk = resk + wgk(2*j)*fc;
    resg = resg + wg(j)*fc; 
  end
  t = wid*xgk(11); 
  fc = feval( f, cen - t ) + feval( f, cen + t );
  resk = wid*( resk + wgk(11)*fc );
  err = abs( resk - wid*resg  );
  return
%
%     end krnrdt
%
