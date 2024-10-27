function output = g_func(x)

k = 6; % case 1

x1 = x(:,1);
x2 = x(:,2);

g1 = 3 + 0.1*(x1-x2).^2 - (x1+x2)/sqrt(2);
g2 = 3 + 0.1*(x1-x2).^2 + (x1+x2)/sqrt(2);
g3 = (x1-x2) + k/sqrt(2);
g4 = (x2-x1) + k/sqrt(2);


output = min([g1,g2,g3,g4],[],2);

end