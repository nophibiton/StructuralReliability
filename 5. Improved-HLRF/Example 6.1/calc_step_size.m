function lambda = calc_step_size(u,Gu,dGu,d, lsffunc,lsfparams,distparams,probdata,Lo)

% assign c > norm(u)/norm(dG(u) )
c      = (norm(u)/norm(dGu))*2 + 10;

% calculate merit at point u
merit = 0.5*(norm(u))^2 + c*abs(Gu);

lambda = 1;
while true
    unew      = u + lambda*d;
    xnew      = trans_u2x(unew,probdata, Lo);

    isGrad    = false;
    Gunew     = eval_lsf(xnew, isGrad, lsffunc,lsfparams,distparams);
    merit_new = 0.5*(norm(unew))^2 + c*abs(Gunew);

    if merit > merit_new
        break;
    end
    lambda = lambda*0.5;
    
end


end