function beta_sys = scm(beta, rho, sys_type, varargin)
%  Performs sequential compounding method for system reliability analysis
%  (Kang et al., 2010)
% 
%  Arguments:
%       beta -- vector of reliability indices of component events, vector
%           size Ncomp x 1
%       rho -- correlation coefficient matrix of all the component events,
%           matrix size Ncomp x Ncomp
%       system type -- string that defines the type of the system which can
%           be: series, parallel, general 
%       additional arguments:
%           system definition - array that contains a vector for the system
%           definition and a string whether a link-set or cut-set
%
%           Example:
%               Ʌ : Intersection, V : Union
%               Esys = E1 Ʌ E2 Ʌ (E3 V E4) Ʌ E5 --> sys_def = {[0 1 0 2 0 3 4 0 5 0],'link'};
%               Esys = (E2 Ʌ E3) V (E1 Ʌ E4) --> sys_def = {[0 2 3 0 1 4 0],'cut'};
%
%
% Returns:
%       beta_sys -- reliability index of the system
%

switch lower(sys_type)

    case 'parallel'

        for m = 1: length(beta)-2
            [beta,rho] = compound_events(beta, rho,'parallel');
        end
        beta_sys = -norminv(mvncdf(-beta,zeros(1,2),rho));

    case 'series'

        for m = 1: length(beta)-2
            [beta,rho] = compound_events(beta, rho,'series');
        end
        beta_sys = norminv(mvncdf(beta,zeros(1,2),rho) );     

    case 'general'
        sys_def = varargin{1};
        sys_ind = sys_def{1};

        E_sys = sys_ind;
        E_sys(E_sys==0) = [];

        psd_R = ones(length(E_sys),length(E_sys) );
        psd_beta = zeros(1,length(E_sys));

        for ik=1:length(E_sys)
            psd_beta(ik) = beta(E_sys(ik) ); 
            for jk=ik:length(E_sys)
                psd_R(ik,jk) = rho(E_sys(ik),E_sys(jk) ); 
                psd_R(jk,ik) = rho(E_sys(ik),E_sys(jk) );
            end
        end
        
        genset_type = sys_def{2};

        switch lower(genset_type)
            case 'link'
                sys_z=find(sys_def{1}==0);
                int_1=sys_z-[0 sys_z(1:end-1)];
                size_sets=int_1(int_1>1)-1;
                Nsets = length(size_sets);

                for i=1:Nsets
                    min_set = sys_def{1}(sys_z(i)+1:sys_z(i)+size_sets(i) );
                    if length(min_set)>1
                        for m = 1:length(min_set)-1
                            [psd_beta,psd_R,E_sys] = compound_events(psd_beta, ...
                                psd_R,'series',[i,i+1], E_sys, rho );
                        end
                    end
                end

                for m = 1: Nsets-2
                    [psd_beta,psd_R,E_sys] = compound_events(psd_beta, psd_R,'parallel',[1,2],E_sys, rho);
                end
                beta_sys = -norminv(mvncdf(-psd_beta,zeros(1,2),psd_R) );  
            
            case 'cut'
          
        end

        
end


