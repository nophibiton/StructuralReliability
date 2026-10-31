   function p = bvnl( dh, dk, r )
%
%  A function for computing bivariate normal probabilities.
%  bvnl calculates the probability that x < dh and y < dk. 
%    parameters  
%      dh 1st upper integration limit
%      dk 2nd upper integration limit
%      r   correlation coefficient
%  Example use:
%   >> bvnl( 3, 1, .35 )
%
%   Author
%       Alan Genz
%       Department of Mathematics
%       Washington State University
%       Pullman, Wa 99164-3113
%       Email : alangenz@wsu.edu
%   This function is based on the method described by 
%        Drezner, Z and G.O. Wesolowsky, (1989),
%        On the computation of the bivariate normal inegral,
%        Journal of Statist. Comput. Simul. 35, pp. 101-107,
%    with major modifications for double precision, for |r| close to 1,
%    and for matlab by Alan Genz - last modifications 7/98.
%
      p = bvnu( -dh, -dk, r );
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
%
%   end bvnu
%
function p = phid(z)
  p = erfc( -z/sqrt(2) )/2;
%
% end phid
%


