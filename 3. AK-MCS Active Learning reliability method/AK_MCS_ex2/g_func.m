function output = g_func(x)


temp = 0;
for i=1:2
    temp = temp + x(:,i).^2 - 5*cos(2*pi().*x(:,i) );
end

output = 10 - temp;

end