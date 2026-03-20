function [gx,dgx] = eval_lsf(x, isGrad, lsffunc,lsfparams,distparams)


func = lsffunc;

% evaluate limit state function at x
gx   = func(x,lsfparams);

% calculate gradient using FFD
dgx = zeros(size(x));
if isGrad
    for i = 1:length(x)
        delta  = distparams(i,2)/200;
        dx     = x;
        dx(i)  = dx(i) + delta;
        gdx    = func(dx,lsfparams);
        dgx(i) = (gdx-gx)/delta;
    end
end

end