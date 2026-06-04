function [g_val,dg_val] = g_fun(x,params)

    R = x(1);
    S = x(2);
    
    g_val = R-S;


    %%% derivative of g(X,params) if available
    dg_val = [];

end