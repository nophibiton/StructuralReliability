% 
% 14.3.1.1 Example: Reliability Analysis of a Column by FORM
%
% Der Kiureghian, A. (2005). First- and Second-Order Reliability Methods. 
% In Engineering Design Reliability Handbook. CRC press.
%
clc,clear; format compact; format shortG;

addpath(genpath(strcat(pwd,'\ferumcore'  )) );

% base units
kN  = 1;
m   = 1;
kPa = kN/m^2;
MPa = 1000*kPa;

% Define Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params       = struct;
params.theta = 2;
params.A     = 0.190 *m^2;
params.S1    = 0.030 *m^3;
params.S2    = 0.015 *m^3;

% Definition of PDF
probdata.marg(1,:) = [ 1  250*kN*m  0.3*250*kN*m  250*kN*m 0 0 0 0 0]; % m1
probdata.marg(2,:) = [ 1  125*kN*m  0.3*125*kN*m  125*kN*m 0 0 0 0 0]; % m2
probdata.marg(3,:) = [ 15 2500*kN   0.2*2500*kN   2500*kN  0 0 0 0 0]; % p
probdata.marg(4,:) = [ 16 40*MPa    0.1*40*MPa    40*MPa   0 0 0 0 0]; % y

% Define the correlation matrix
probdata.correlation   = [1.0,0.5,0.3,0.0; ...
                          0.5,1.0,0.3,0.0; ...
                          0.3,0.3,1.0,0.0; ...
                          0.0,0.0,0.0,1.0];
probdata.parameter     = distribution_parameter(probdata.marg);

% Define transformation using Nataf model
Ro  = mod_corr(probdata, probdata.correlation );
Lo  = (chol(Ro))'; % Cholesky decomposition
u2x = @(u) u_to_x(u', probdata,Lo);

% Define limit state function in G(T(u))
G_FUN = @(u) g_fun(u2x(u),params); 

N        = 1000;
p0       = 0.1;
max_it   = 50;

while true

Nc       = N*p0;             % number of markov chains
g_vals   = zeros(1,N);


% (1) generate iid samples from standard normal space phi_n
u_j = randn(size(probdata.marg,1), N)';

% (2) Order the samples increasing order of  magnitude of their limit-state values
%     Find c1 as the p0-percentile of the samples. Fj = {G(T(u) ) <= cj}
for k = 1:N
   g_vals(k) = G_FUN(u_j(k,:));
end

% (3) Iterations
j    = 1;

    
% (4) Repeat while cj > 0
while true
        % determine p0 percentile value
        c(j) = prctile(g_vals, p0*100);
        
        % (4a) Generate N samples using MCMC
        % select seeds for MCMC sampling included in {uj, G(uj)<= cj}
        [g_vals_sorted, idx]   = sort(g_vals);
        u_j_sorted             = u_j(idx,:); 
        nF                     = sum(g_vals <= max(c(j), 0));
        u_j_rnd_seeds          = u_j_sorted(randperm(nF),:); % randomized order
        
        % assign conditional probability to the level
        if c(j) <= 0
            c(j)    = 0;
            prob(j) = nF/N;
        else
            prob(j) = p0;      
        end

   % compute coefficient of variation
   if j == 1
      delta(j) = sqrt(((1-p0)/(N*p0)));   % cov for p(1): MCS (Ref. 2 Eq. 8)
   else
      I_Fj     = reshape(g_vals <= c(j),(1/p0),Nc);         % indicator function for the failure samples
      p_j      = (1/N)*sum(I_Fj(:));                   % ~=p0, sample conditional probability
      gamma    = corr_factor(I_Fj,p_j,(1/p0),Nc);          % corr factor (Ref. 2 Eq. 10)
      delta(j) = sqrt( ((1-p_j)/(N*p_j))*(1+gamma) );  % coeff of variation(Ref. 2 Eq. 9)
   end

        
        % initialize parameters of the values to be sampled
        Ns      = size(u_j_rnd_seeds, 1); % number of seeds
        n       = size(u_j_rnd_seeds, 2); % dimension of random vector
        u_j_gen = zeros(N,n);             % generated samples
        g_j_gen = zeros(1,N);             % LSF evaluations of G(T(uj))
        
        % define the number of chains and spread/parameter of MCMC
        Nchain              = ones(Ns,1)*floor(N/Ns);
        Nchain(1:mod(N,Ns)) = Nchain(1:mod(N,Ns))+1;  % number of chains per seed
        sigma               = std(u_j_rnd_seeds,0,1); % spread of MCMC per seed
        
        % begin sampling per seed
        for k = 1:Ns
            idx            = sum(Nchain(1:k-1))+1;
            u_j_gen(idx,:) = u_j_rnd_seeds(k,:);
            g_j_gen(idx)   = G_FUN(u_j_gen(idx,:));
            for t = 1:Nchain(k)-1 % iterate number of chains per seed
                % M-H component by component
                u0 = u_j_gen(idx+t-1,:);
                for i = 1:n
                    % 1a. generate a pre-candidate xi by sampling from q(|u0i)
                    xi(i) = u0(i)+(2*rand()-1)*sigma(i);
                    % 1b. accept or reject xi
                    ai = min(1,normpdf(xi(i))/normpdf(u0(i)) );
                    if rand() < ai
                        v(i) = xi(i);
                    else
                        v(i) = u0(i);
                    end
                end
                % 2. accept or reject v
                temp_geval = G_FUN(v);
                if temp_geval <= c(j)
                    % accept new state
                    u_j_gen(idx+t,:) = v;
                    g_j_gen(idx+t  ) = temp_geval;
                else
                    % reject and use previous state
                    u_j_gen(idx+t,:) = u_j_gen(idx+t-1,:);
                    g_j_gen(idx+t  ) = g_j_gen(idx+t-1);
                end
            end
        end
        
        u_j    = u_j_gen;
        g_vals = g_j_gen;
        
        j = j+1;
        

        if c(j-1) <= 0 || j-1 == max_it
           break;
        end

end

Pf = prod(prob)

% coeficient of variation estimate
delta_SuS = sqrt(sum(delta.^2))   % (Ref. 2 Eq. 12)

if delta_SuS <= 0.05 
    break;
else
    N = N + 10000;
end

end

% Pfs(it) = Pf;
% end
% 
% hist(Pfs,50)
% mean(Pfs)

rmpath(genpath(strcat(pwd,'\ferumcore'))   );