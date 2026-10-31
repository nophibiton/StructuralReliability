%
%  This file contains functions tvnls (trivariate normal), bvnl 
%   (bivariate normal), phid (univariate normal), plus support functions.
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
   function tvn = tvnls( h, r )
%    
%     A function for computing trivariate normal probabilities.
%     It calculates the probability that x(i) < h(i), for i = 1, 2, 3.     
%    h     real array of three upper limits for probability distribution 
%    r     real array of three correlation coefficients, r should 
%          contain the lower left portion of the correlation matrix R; 
%          r should contain the values r21, r31, r32 in that order.
% 
%     This function uses algorithms developed from the ideas 
%     described in the papers:
%       R.L. Plackett, Biometrika 41(1954), pp. 351-360.
%       Z. Drezner, Math. Comp. 62(1994), pp. 289-294.
%     with Gauss-Legendre rule integration from (0,0,r32) to R. 
%     Absolute accuracy for typical calculations is 14-15 digits.
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
      epst = 1e-14;
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
      tvn = 0;
      if ( abs(h1) + abs(h2) + abs(h3) < epst )
         tvn = ( 1 + ( asin(r12) + asin(r13) + asin(r23) )/pt )/8;
      elseif ( abs(r12) + abs(r13) < epst )
         tvn = phid(h1)*bvnl( h2, h3, r23 );
      elseif ( abs(r13) + abs(r23) < epst ) 
         tvn = phid(h3)*bvnl( h1, h2, r12 );
      elseif ( abs(r12) + abs(r23) < epst ) 
         tvn = phid(h2)*bvnl( h1, h3, r13 );
      elseif ( 1 - r23 < epst ) 
         tvn = bvnl( h1, min( [h2 h3] ), r12 );
      elseif ( r23 + 1 < epst ) 
         if  ( h2 > -h3 ) 
             tvn = bvnl( h1, h2, r12 ) - bvnl( h1, -h3, r12 );
         end
      else
%
%       Compute singular TVN value
%
        tvn = bvnl( h2, h3, r23 )*phid(h1);
%
%       Use Gauss-Legendre numerical integration to compute probability
%
        wg( 1: 3) = [0.1279381953467518 0.1258374563468280 0.1216704729278031];
        wg( 4: 6) = [0.1155056680537265 0.1074442701159659 .09761865210411358];
        wg( 7: 9) = [.08619016153195296 .07334648141108081 .05929858491543594];
        wg(10:12) = [.04427743881742087 .02853138862893389 .01234122979998693];
        xg( 1: 3) = [.06405689286260559 0.1911188674736164 0.3150426796961635];
        xg( 4: 6) = [0.4337935076260450 0.5454214713888396 0.6480936519369754];
        xg( 7: 9) = [0.7401241915785546 0.8200019859739028 0.8864155270044012];
        xg(10:12) = [0.9382745520027329 0.9747285559713096 0.9951872199970214];
        rua = asin( r12 ); rub = asin( r13 ); res = 0;
        for j = 1 : 12
          fc = 0;
          [ r12t rr2 ] = sincs( rua*( 1 - xg(j) )/2 ); 
          [ r13t rr3 ] = sincs( rub*( 1 - xg(j) )/2 );
          if ( abs(rua) > 0 ) 
             fc = fc + rua*pntgnd( h1, h2, h3, r13t, r23, r12t, rr2 );
          end
          if ( abs(rub) > 0 )
             fc = fc + rub*pntgnd( h1, h3, h2, r12t, r23, r13t, rr3 );
          end
          [ r12t rr2 ] = sincs( rua*( 1 + xg(j) )/2 ); 
          [ r13t rr3 ] = sincs( rub*( 1 + xg(j) )/2 );
          if ( abs(rua) > 0 ) 
             fc = fc + rua*pntgnd( h1, h2, h3, r13t, r23, r12t, rr2 );
          end
          if ( abs(rub) > 0 )
             fc = fc + rub*pntgnd( h1, h3, h2, r12t, r23, r13t, rr3 );
          end
          res = res + wg(j)*fc; 
        end
        tvn = tvn + res/( 4*pi );
      end
      tvn = max( 0, min( tvn, 1 ) ) ;
   return
%
% end tvnls
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
   function f = pntgnd( ba, bb, bc, ra, rb, r, rr )
%
%     Computes Plackett formula integrand;
%
      f = 0;
      dt = rr*( rr - ( ra - rb )^2 - 2*ra*rb*( 1 - r ) );
      if ( dt > 0 ) 
         bt = ( bc*rr + ba*( r*rb - ra ) + bb*( r*ra -rb ) )/sqrt(dt) ;
         ft = ( ba - r*bb )^2/rr + bb*bb;
         if ( bt > -10 & ft < 100 ) 
            f = exp( -ft/2 );
            if ( bt < 10 ),  f = f*phid(bt); end
         end
      end
   return
%
% end pntgnd
%
   function p = phid(z)
      p = erfc( -z/sqrt(2) )/2;
   return
%
% end phid
%
   function p = bvnl( dh, dk, r )
%
%  A function for computing bivariate normal probabilities.
%  bvnl calculates the probability that x < dh and y < dk. 
%    parameters  
%      dh 1st upper integration limit
%      dk 2nd upper integration limit
%      r   correlation coefficient
%
%   Author
%       Alan Genz
%       Department of Mathematics
%       Washington State University
%       Pullman, Wa 99164-3113
%       Email : alangenz@wsu.edu
%   This function is based on the method described by 
%        Drezner, Z and G.O. Wesolowsky, (1989),
%         On the computation of the bivariate normal inegral,
%          Journal of Statist. Comput. Simul. 35, pp. 101-107,
%    with major modifications for double precision, for |r| close to 1,
%    and for matlab by Alan Genz - last modifications 7/98.
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






