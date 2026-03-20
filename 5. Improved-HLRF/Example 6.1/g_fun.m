function [g_val,dg_val] = g_fun(x,params)

    X1 = x(1);
    X2 = x(2);
    
    g_val = X1^2 - 2*X2;


    %%% derivative of g(X,params) if available
    dg_val = [];

end