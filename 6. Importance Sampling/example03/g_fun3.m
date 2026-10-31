function [g_val,dg_val] = g_fun3(x,params)

    M1 = x(1);
    M2 = x(2);    
    M3 = x(3);
    M4 = x(4);    
    M5 = x(5);
    H  = x(6);    
    V  = x(7);   

    g_val = M1 + 2*M3 + 2*M4 + M5 - 5*H - 5*V;

    %%% derivative of g(X,params) if available
    dg_val = [];

end