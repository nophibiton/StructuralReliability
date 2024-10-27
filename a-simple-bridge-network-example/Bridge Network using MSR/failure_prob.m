function [Pf,varargout]=failure_prob(beta,r,c_sys,sys_type,sys_def,integ_method,varargin)

m=size(r,2);
if isempty(varargin)~=1;
    SORM_method=varargin{1};
    if strcmpi(SORM_method,'point fitting')==1;
        point_fitting_method=varargin{2};
    end
end
        
switch lower(integ_method)
    case 'direct'
        tol=10^-7; int_lb=-6; int_ub=6;
        switch m
            case 1
                Pf=quad(@(x)integrnd(beta,r,c_sys,sys_type,x),...
                    int_lb,int_ub,tol);
            case 2
                Pf=dblquad(@(x,y)integrnd(beta,r,c_sys,sys_type,x,y),...
                    int_lb,int_ub,int_lb,int_ub,tol);
            case 3
                Pf=triplequad(@(x,y,z)integrnd(beta,r,c_sys,sys_type,...
                    x,y,z),int_lb,int_ub,int_lb,int_ub,int_lb,int_ub,tol);
            case 4
                %%%%%%%%%%%
                Pf = 0;
                X_range = 6; %X_range = 4;
                interval = 0.2; %interval = 0.2;
                for X1=-X_range:interval:X_range
                    for X2=-X_range:interval:X_range
                        for X3=-X_range:interval:X_range
                            for X4=-X_range:interval:X_range
                                % approximate P(Ei|S=s) using DS-model
                                p_s = normcdf((-beta-r(:,1)*X1-r(:,2)*X2-r(:,3)*X3-r(:,4)*X4)./sqrt(1-r(:,1).^2-r(:,2).^2-r(:,3).^2-r(:,4).^2));
                                % expand into MCE events
                                p_s = prob_vector(p_s); 
                                % calculate system reliability using numerical integration
                                Pf = Pf + c_sys'*p_s*normpdf(X1)*normpdf(X2)*normpdf(X3)*normpdf(X4)*interval^4;
                            end
                        end
                    end
                end 
                % added by NIB 2023
                %%%%%%%%%%

        end
    case {'form','sorm'}
        [probdata,analysisopt]=create_input(m);
        formresults=form_msr(c_sys,sys_type,beta,r,probdata,analysisopt);
        if isempty(formresults)==1; 
            Pf=[];
            if strcmpi(integ_method,'SORM')==1;
                varargout={[],[],[]};
            end
            return;
        end
        Pf=formresults.pf1;
        if strcmpi(integ_method,'SORM')==1;
            switch lower(SORM_method);
                case 'curvature fitting'
                    sorm_curfit_results=sorm_curvature_fitting_msr(c_sys,sys_type,beta,r,formresults,probdata);
                    if isempty(sorm_curfit_results)==1; varargout={[],[],[]}; return; end;
                    Pf_Breitung=sorm_curfit_results.pf2_breitung;
                    Pf_imp_Breitung=sorm_curfit_results.pf2_breitung_mod;
                    Pf_Tvedt=sorm_curfit_results.pf2_tvedt_EI;
                case 'point fitting'
                    switch lower(point_fitting_method)
                        case 'secant'
                            sorm_point_fitted_secant_results=sorm_point_fitted_mod_secant_msr(c_sys,sys_type,beta,r,formresults,probdata,analysisopt.ig_max);
                            if isempty(sorm_point_fitted_secant_results)==1; varargout={[],[],[]}; return; end;
                            Pf_Breitung=sorm_point_fitted_secant_results.pf2_breitung;
                            Pf_imp_Breitung=sorm_point_fitted_secant_results.pf2_breitung_mod;
                            Pf_Tvedt=sorm_point_fitted_secant_results.pf2_tvedt_EI;
                        case 'newton'
                            sorm_point_fitted_newton_results=sorm_point_fitted_mod_msr(c_sys,sys_type,beta,r,formresults,probdata,analysisopt.ig_max);
                            if isempty(sorm_point_fitted_newton_results)==1; varargout={[],[],[]}; return; end;
                            Pf_Breitung=sorm_point_fitted_newton_results.pf2_breitung;
                            Pf_imp_Breitung=sorm_point_fitted_newton_results.pf2_breitung_mod;
                            Pf_Tvedt=sorm_point_fitted_newton_results.pf2_tvedt_EI;
                    end
            end
        end
end

if strcmpi(sys_type,'series')==1 || (strcmpi(sys_type,'general')==1 && strcmpi(sys_def{2},'link')==1);
    Pf=1-Pf;
    if strcmpi(integ_method,'SORM')==1;
        Pf_Breitung=1-Pf_Breitung; Pf_imp_Breitung=1-Pf_imp_Breitung; Pf_Tvedt=1-Pf_Tvedt;        
    end
end
if isempty(varargin)~=1;
    varargout={Pf_Breitung,Pf_imp_Breitung,Pf_Tvedt};
end

end
