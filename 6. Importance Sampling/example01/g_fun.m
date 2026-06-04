function [g_val,dg_val] = g_fun(x,params)

    m1 = x(1);
    m2 = x(2);
    p  = x(3);
    y  = x(4);
    
    theta = params.theta;
    A     = params.A;
    S1    = params.S1;
    S2    = params.S2;
    
    
    g_val = 1 - (m1./(S1*y)) - (m2./(S2.*y)) - (p./(A.*y)).^theta;


    %%% derivative of g(X,params) if available
    dg_val = [];

end