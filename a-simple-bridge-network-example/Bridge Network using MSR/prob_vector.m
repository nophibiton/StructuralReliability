function [p]=prob_vector(p_marginals)

for i=1:length(p_marginals);
    if i==1;
        p=[p_marginals(1); (1-p_marginals(1))];
    else
        p=[p*p_marginals(i); p*(1-p_marginals(i))];
    end
end

end