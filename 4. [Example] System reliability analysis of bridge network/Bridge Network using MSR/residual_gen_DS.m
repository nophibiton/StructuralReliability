function norm_R=residual_gen_DS(r,R)

n=length(R);
norm_R=norm(R-(r*r'-diag(diag(r*r'))+eye(n)),'fro');

end