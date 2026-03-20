function Ro = mod_corr(R,probdata)

% unpack
Xdist = probdata.Xdist;
mu    = probdata.Xmu;
std   = probdata.Xstd;

nrv = numel(mu);
Ro  = zeros(nrv,nrv);

for i = 1:nrv
    for j = i+1:nrv
        switch lower(Xdist{i})
            %%%%%%%%%%%%%
            case 'normal'
                switch lower(Xdist{j})
                    case 'normal' 
                        Ro(i,j) = R(i,j);
                    case 'lognormal'
                        delta_j = std(j)/mu(j);
                        C       = delta_j/sqrt(log(1+delta_j^2));
                        Ro(i,j) = C*R(i,j);
                    case 'gumbel' % type I largest
                        C       = 1.031;
                        Ro(i,j) = C*R(i,j);
                    case 'weibull' % type III smallest (epsilon=0)
                        delta_j = std(j)/mu(j); 
                        C       = 1.031 - 0.195*delta_j + 0.328*delta_j^2;
                        Ro(i,j) = C*R(i,j);
                end
            %%%%%%%%%%%%% 
            case 'lognormal'
                switch lower(Xdist{j})
                    case 'normal'
                        delta_i = std(i)/mu(i);
                        C       = delta_i/sqrt(log(1+delta_i^2));
                        Ro(i,j) = C*R(i,j);
                    case 'lognormal'
                        delta_i = std(i)/mu(i);
                        delta_j = std(j)/mu(j);
                        C       = log(1+R(i,j)*delta_i*delta_j)/...
                                  (R(i,j)*sqrt(log(1+delta_i^2)*... 
                                  log(1+delta_j^2) ) );
                        Ro(i,j) = C*R(i,j);
                    case 'gumbel' % type I largest
                        delta_i = std(i)/mu(i);
                        C       =  1.029 + 0.001*R(i,j) + 0.014*delta_i +...
                                   0.004*R(i,j)^2 + 0.233*delta_i^2 - ...
                                   0.197*R(i,j)*delta_i ;
                        Ro(i,j) = C*R(i,j);
                    case 'weibull' % type III smallest (epsilon=0)
                        delta_i = std(i)/mu(i); 
                        delta_j = std(j)/mu(j);
                        C       = 1.031 + 0.052*R(i,j) + 0.011*delta_j -...
                                  0.210*delta_i + 0.002*R(i,j)^2 + ...
                                  0.220*delta_j^2 + 0.350*delta_i^2 + ...
                                  0.009*delta_i*delta_j + 0.005*R(i,j)*...
                                  delta_j - 0.174*R(i,j)*delta_i;
                        Ro(i,j) = C*R(i,j);
                end
            %%%%%%%%%%%%% 
            case 'gumbel' % type I largest
                switch lower(Xdist{j})
                    case 'normal'
                        C       = 1.031;
                        Ro(i,j) = C*R(i,j);
                    case 'lognormal'
                        delta_j = std(j)/mu(j);
                        C       = 1.029 + 0.001*R(i,j) + 0.014*delta_j +...
                                  0.004*R(i,j)^2 + 0.233*delta_j^2 - ...
                                  0.197*R(i,j)*delta_j;
                        Ro(i,j) = C*R(i,j);
                    case 'gumbel' % type I largest
                        C       = 1.064 - 0.069*R(i,j) + 0.005*R(i,j)^2 ;
                        Ro(i,j) = C*R(i,j);
                    case 'weibull' % type III smallest (epsilon=0)
                        delta_j = std(j)/mu(j);
                        C       = 1.064 + 0.065*R(i,j) - 0.210*delta_j+...
                                  0.003*R(i,j)^2 + 0.356*delta_j^2 - ...
                                  0.211*R(i,j)*delta_j;          
                        Ro(i,j) = C*R(i,j);
                end

            %%%%%%%%%%%%% 
            case 'weibull' % type III smallest (epsilon=0)
                switch lower(Xdist{j})
                    case 'normal'
                        delta_i = std(i)/mu(i);
                        C       = 1.031 - 0.195*delta_i + 0.328*delta_i^2;
                        Ro(i,j) = C*R(i,j);
                    case 'lognormal'
                        delta_j = std(j)/mu(j);
                        delta_i = std(i)/mu(i);
                        C       = 1.031 + 0.052*R(i,j) + 0.011*delta_i - ...
                                  0.210*delta_j + 0.002*R(i,j)^2 + ...
                                  0.220*delta_i^2 + 0.350*delta_j^2 + ...
                                  0.009*delta_i*delta_j + 0.005*R(i,j)*...
                                  delta_i - 0.174*R(i,j)*delta_j;
                        Ro(i,j) = C*R(i,j);
                    case 'gumbel' % type I largest
                        delta_i = std(i)/mu(i);
                        C       =  1.064 + 0.065*R(i,j) - 0.210*delta_i  + ...
                                   0.003*R(i,j)^2 + 0.356*delta_i^2 - ...
                                   0.211*R(i,j)*delta_i;
                        Ro(i,j) = C*R(i,j);                        
                    case 'weibull' % type III smallest (epsilon=0)
                        delta_i = std(i)/mu(i);
                        delta_j = std(j)/mu(j);
                        C       = 1.063 - 0.004*R(i,j) - 0.200*(delta_i+...
                                  delta_j) - 0.001*R(i,j)^2 + 0.337*...
                                  (delta_i^2+delta_j^2) + 0.007*R(i,j)*...
                                  (delta_i+delta_j) - 0.007*delta_j*delta_i;          
                        Ro(i,j) = C*R(i,j);
                end

                             

        end


    end
end

Ro = Ro + Ro.' + eye (nrv);


end