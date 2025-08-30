function [beta,rho, E_sys] = compound_events(beta, rho, type, varargin)
%
%  
%
%
%

if isempty(varargin) 
    pair=[1,2]; 
    E_sys = 1:1:length(beta);
    rho_orig = rho;
else 
    pair = varargin{1}; 
    E_sys = varargin{2};
    rho_orig = varargin{3};
end



switch lower(type)

    case 'series' %% SERIES EVENTS
            c1 = pair(1);
            c2 = pair(2);
            new_combined_ind = c1;
            
            % STEP 1: calculate new beta of compounded event
            beta_1 = beta(c1); beta_2 = beta(c2);

            rho_12 = rho(c1,c2);

            R1 = [1,rho_12;rho_12,1];
            beta_1or2 = norminv(mvncdf([beta_1,beta_2],zeros(1,2),R1));

            
            % replace new betas
            new_beta = beta;
            new_beta(c1) = [];
            new_beta(c1) = beta_1or2;
            

            % STEP 2: calculate new correlation coefficients
            new_rho = rho;
            new_rho(:,c2) = [];  new_rho(c2,:) = [];
            new_rho(:,c1) = 1;   new_rho(c1,:) = 1; 
            new_E_sys = E_sys; 
            new_E_sys(c2) = [];
            new_E_sys(c1) = -1;
            
            % get index of remaining components after compounding
            index_vect = 1:length(beta);
            index_vect(index_vect==c1) = [];
            index_vect(index_vect==c2) = [];
            
            for k = index_vect


                beta_k = beta(k );
                rho_1k = rho(c1,k );
                rho_2k = rho(c2,k );
                
                p1 = [beta_1,beta_2,beta_k,beta_1or2 ];
                p2 = [rho_12,rho_1k,rho_2k ];
                p3 = 0;

                if E_sys(k) == E_sys(c1) || E_sys(k) == E_sys(c2), p3 = 1; end 


                options = optimoptions('fmincon','Display','none', ...
                    'StepTolerance',1e-9,'OptimalityTolerance',1e-9, ...
                    'ConstraintTolerance',1e-9);
                x = fmincon(@(x)obj_fun_union(x,p1,p2,p3),0, ...
                    [],[],[],[],-1,1,[],options );

                if k<new_combined_ind
                    new_rho(k,new_combined_ind) = x; new_rho(new_combined_ind,k) = x;
                else
                    new_rho(k-1,new_combined_ind) = x; new_rho(new_combined_ind,k-1) = x;
                end

            end
        
            % Update values of beta and rho
            beta = new_beta;  
            rho = new_rho;
            E_sys = new_E_sys;


    case 'parallel' %% PARALLEL EVENTS
            c1 = pair(1);
            c2 = pair(2);
            new_combined_ind = c1;
            
            % STEP 1: calculate new beta of compounded event
            beta_1 = beta(c1);
            beta_2 = beta(c2);
            rho_12 = rho(c1,c2);
            R1 = [1,rho_12;rho_12,1];
            beta_1and2 = -norminv(mvncdf([-beta_1,-beta_2],zeros(1,2),R1));

            % replace new betas
            new_beta = beta;
            new_beta(c1) = [];
            new_beta(c1) = beta_1and2;
            
            % STEP 2: calculate new correlation coefficients
            new_rho = rho;
            new_rho(:,c2) = [];  new_rho(c2,:) = [];
            new_rho(:,c1) = 1;   new_rho(c1,:) = 1;
            new_E_sys = E_sys; 
            new_E_sys(c2) = [];
            new_E_sys(c1) = -1;
            
            % get index of remaining components after compounding
            index_vect = 1:length(beta);
            index_vect(index_vect==c1) = [];
            index_vect(index_vect==c2) = [];
            
            for k = index_vect

                beta_k = beta(k);
                rho_1k = rho(c1,k);
                rho_2k = rho(c2,k);

                p1 = [beta_1,beta_2,beta_k,beta_1and2];
                p2 = [rho_12,rho_1k,rho_2k];
                p3 = 0;

                if E_sys(k) == E_sys(c1) || E_sys(k) == E_sys(c2), p3 = 1; end
%                 p3
                options = optimoptions('fmincon','Display','none', ...
                    'StepTolerance',1e-9,'OptimalityTolerance',1e-9, ...
                    'ConstraintTolerance',1e-9);
                x = fmincon(@(x)obj_fun_intersect(x,p1,p2,p3),0.1, ...
                    [],[],[],[],-1,1,[],options);

                if k<new_combined_ind
                    new_rho(k,new_combined_ind) = x; 
                    new_rho(new_combined_ind,k) = x;
                else
                    new_rho(k-1,new_combined_ind) = x; 
                    new_rho(new_combined_ind,k-1) = x;
                end
        
            end

            % Update values of beta and rho
            beta = new_beta ;
            rho = new_rho;
            E_sys = new_E_sys;

end % end of switch


