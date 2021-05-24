function [lambda] = cal_lambda(N, alpha, w, p)
    % w is the price of each input factor
    % p is the price of rice
    lambda=1.0;
    for i=1:N
        lambda=lambda*(alpha(i)/(w(i)/p))^(alpha(i));
    end
    
end

