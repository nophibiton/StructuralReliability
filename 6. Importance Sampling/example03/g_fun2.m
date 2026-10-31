function [g_val,dg_val] = g_fun2(x,params)

    M1 = x(1);
    M2 = x(2);    
    M3 = x(3);
    M4 = x(4);    
    M5 = x(5);
    H  = x(6);    
    V  = x(7);   

    g_val = M2 + 2*M3 + M4 - 5*V;

    %%% derivative of g(X,params) if available
    dg_val = [];

end