function std_dist = calculate_std(sup,prob,n)
    % n is the number of grids in discretization
    % sup is the support of the distribution
    % prob is the probability of each point in the support
    std_dist=0.0;
    avg=dot(sup,prob);
    for i=1:n
        std_dist=std_dist+(sup(i)-avg)^2*prob(i);
    end
    std_dist=std_dist^0.5;


end

