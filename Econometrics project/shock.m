function [area_sup, cec_sup, organic_sup, area_prob, cec_prob, organic_prob,...
        eps0_sup, eps1_sup, eps2_sup, eps3_sup, eps0_prob, eps1_prob, eps2_prob, eps3_prob,...
        mean_initial, mean_eps0, mean_eps1, mean_eps2, mean_eps3,...
        mean_exp_eps0_combination, mean_exp_eps1_combination,...
        mean_exp_eps2_combination, mean_exp_eps3_combination] = shock(n)
	% this function constructs the distribution of the shocks
    
    [area_sup, area_prob]=ini_normal(1.12, 1.92, n);
    calculate_std(area_sup,area_prob,n);
    [cec_sup, cec_prob]=ini_normal(2.32, 1.20, n);
    calculate_std(cec_sup,cec_prob,n);
    [organic_sup, organic_prob]=ini_normal(0.52, 0.36, n);
    calculate_std(organic_sup,organic_prob,n);

    [eps0_sup(1,:), eps0_prob(1,:)]=ini_normal(49.94, 25.48, n); % period 0 rainfall
    calculate_std(eps0_sup(1,:),eps0_prob(1,:),n);
    [eps0_sup(2,:), eps0_prob(2,:)]=ini_normal(33.85, 0.49, n); % period 0 temp
    calculate_std(eps0_sup(2,:),eps0_prob(2,:),n);

    [eps1_sup(1,:), eps1_prob(1,:)]=ini_normal(25.55, 13.28, n); % period 1 rainfall
    calculate_std(eps1_sup(1,:),eps1_prob(1,:),n);
    [eps1_sup(2,:), eps1_prob(2,:)]=ini_normal(32.95, 0.91, n); % period 1 temp
    calculate_std(eps1_sup(2,:),eps1_prob(2,:),n);
    
    [eps2_sup(1,:), eps2_prob(1,:)]=ini_normal(29.24, 16.24, n); % period 2 rainfall
    calculate_std(eps2_sup(1,:),eps2_prob(1,:),n);
    [eps2_sup(2,:), eps2_prob(2,:)]=ini_normal(32.19, 0.67, n); % period 2 temp
    calculate_std(eps2_sup(2,:),eps2_prob(2,:),n);
    
    [eps3_sup(1,:), eps3_prob(1,:)]=ini_normal(2.92, 4.96, n); % period 3 rainfall
    calculate_std(eps3_sup(1,:),eps3_prob(1,:),n);
    [eps3_sup(2,:), eps3_prob(2,:)]=ini_normal(31.09, 1.32, n); % period 3 temp
    calculate_std(eps3_sup(2,:),eps3_prob(2,:),n);
    
    mean_initial(1:3)=0.0;
    mean_initial(1)=dot(area_sup', area_prob');
    mean_initial(2)=dot(cec_sup', cec_prob');
    mean_initial(3)=dot(organic_sup', organic_prob');
    
    mean_eps0=dot(eps0_sup', eps0_prob');
    mean_eps1=dot(eps1_sup', eps1_prob');
    mean_eps2=dot(eps2_sup', eps2_prob');
    mean_eps3=dot(eps3_sup', eps3_prob');
    
    mean_exp_eps0_combination=mean_exp_eps_combination(n, eps0_sup, eps0_prob, mean_eps0);
    mean_exp_eps1_combination=mean_exp_eps_combination(n, eps1_sup, eps1_prob, mean_eps1);
    mean_exp_eps2_combination=mean_exp_eps_combination(n, eps2_sup, eps2_prob, mean_eps2);
    mean_exp_eps3_combination=mean_exp_eps_combination(n, eps3_sup, eps3_prob, mean_eps3);
    
end

