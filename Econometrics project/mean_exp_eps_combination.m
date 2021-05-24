%%%% For example, this files calculate E(exp(eps3)), eps3 is constructed
%%%% from rainfall3 and temperature3, as
%%%% 0.5*(2-rainfall/mean(rainfall)-temp/mean(temp))
function [v] = mean_exp_eps_combination(n, eps_sup, eps_prob, mean_eps)
    
    % this function calculates epsilon with realized rain and temperature
    v=0;
    for i_rain=1:n
        for i_temp=1:n
            v=v+exp((2-(abs(eps_sup(1,i_rain)-mean_eps(1))/mean_eps(1)+...
                abs(eps_sup(2,i_temp)-mean_eps(2))/mean_eps(2)))*0.5)...
                *eps_prob(1,i_rain)*eps_prob(2,i_temp);
        end
    end
    
end

