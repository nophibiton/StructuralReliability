function [C]=event_matrix(n_events)
%
% This function creates a recursive formulation of making the 
% event matrix of all MECE events given the number of component 
% events Ncomp.
%

for i=1:n_events;
    if i==1;
        C=[1; 0];
    else
        C=[C ones(length(C),1); C zeros(length(C),1)];
    end
end

end