clc, clear; format compact;

n_events = 7; % number of component events

% defined the correlation coefficient
r = [0.90, 0.96, 0.91, 0.95, 0.92, 0.94, 0.93];
R = r'*r; for i=1:size(R,1), R(i,i)=1; end

% calculate the unicomponent probability in Eq. 11
ui = (100/sqrt(3)-200)/40;
P_Ei = mvncdf (ui,0);

% define the event matrix by zeros an ones 
[C] = event_matrix(n_events); % procedure is defined by Kang 2012
n_mece = 2^n_events;

% Define the problem data
c_sys = [ones(n_mece-1,1);0]; % define the event vector for a series event

% Solve the problem using CVX for the lower bound
cvx_begin
    variable p(n_mece)       
    minimize(c_sys' * p)    
    subject to

        sum(p) == 1      
        for i=1:n_mece
            p(i) >= 0 
        end

        for i=1:n_events
            C(:,i)'*p==P_Ei
        end
 
cvx_end

LowerBound = cvx_optval

% Solve the problem using CVX for the upper bound
cvx_begin
    variable p(n_mece)       
    maximize(c_sys' * p)    
    subject to

        sum(p) == 1      
        for i=1:n_mece
            p(i) >= 0 
        end

        for i=1:n_events
            C(:,i)'*p==P_Ei
        end
 
cvx_end

UpperBound = cvx_optval
