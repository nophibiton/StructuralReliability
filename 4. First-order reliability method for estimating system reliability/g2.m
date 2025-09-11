function val = g2(x)

D = 610;
t = 7.16;
l = 50;

% sigma_u = x(5,:);
% p = x(6,:);
% d = x(7,:);
% e = x(8,:);

sigma_u = x(1,:);
p = x(2,:);
d = x(3,:);
e = x(4,:);


q = l^2/(D*t);
if q > 50   
    M = 3.3 + 0.032*q;
else
    M = sqrt(1 + 0.6275*q - 0.003375*q^2 );
end

pb = e*(1.8*t*sigma_u/D) * ((1-d/t)/(1-d/(M*t)) );

% limit state function
val = pb - p;

end