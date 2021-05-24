function [x1] = solve_x1(g1,gamma, theta, alpha1, alpha2, alpha3, A, B, ratio, ...
                        eps0_sup, eps0_prob, mean_eps0, ...
                        mean_exp_eps1_combination, y0_real, eps0_real, price)
    % this function solves stage 1 input values with given structural parameter values
    [n, obs, ~, ~, ~, M, N, ~, ~, ~, ~, ~,...
    ~, ~,...
    ~, ~, ~,~]=parameter();

    lambda1=cal_lambda(N(1), alpha1, price(1:N(1)), price((N(1)+N(2)+N(3))+1));
    lambda2=cal_lambda(N(2), alpha2, price((N(1)+1):(N(1)+N(2))), price((N(1)+N(2)+N(3))+1));
    lambda3=cal_lambda(N(3), alpha3, price((N(1)+N(2)+1):(N(1)+N(2)+N(3))), price((N(1)+N(2)+N(3))+1));
    
    mean_w1=price(1:N(1));
    mean_w2=price((N(1)+1):(N(1)+N(2)));
    mean_w3=price((N(1)+N(2)+1):(N(1)+N(2)+N(3)));
    mean_p=price((N(1)+N(2)+N(3))+1);
    
    x1(1:N(1))=0.0;
    
    deviation=[abs(eps0_real(1)-mean_eps0(1))/mean_eps0(1),...
                abs(eps0_real(2)-mean_eps0(2))/mean_eps0(2)];
    
    EP2=0.0; %% Don't forget to time theta2^(1/gamma2)
    value3=0.0; %% store intermediate value for Q3^((gamma3-1)/gamma3)
    value2=0.0; %% stode intermediate value for Q2^((gamma2-1)/gamma2)
    eps3_combination=1-ratio(4,:)*deviation';
    eps2_combination=1-ratio(3,:)*deviation';
    eps1_combination=1-ratio(2,:)*deviation';
    eps0_combination=1-ratio(1,:)*deviation';
    
%     value3=((A(3)*exp(eps3_combination))^(gamma(3)/(gamma(3)-1))-...
%            ((1-theta(3))*(B(3)*lambda3)^gamma(3))^(1/(1-gamma(3))))^((gamma(3)-1)/gamma(3));
%        
%     value2=((A(2)*exp(eps2_combination)*theta(3)^(1/gamma(3))*value3)^(gamma(2)/(gamma(2)-1))-...
%            ((1-theta(2))*(B(2)*lambda2)^gamma(2))^(1/(1-gamma(2))))^((gamma(2)-1)/gamma(2));
%     EP2=value2*theta(2)^(1/gamma(2));
%     g1=(1-theta(1))^(1/(gamma(1)-1))*...
%        (A(1)*B(1)*lambda1*EP2*exp(eps1_combination))^(gamma(1)/(gamma(1)-1))-1;

    left=(theta(1)/(1-theta(1))/g1)^(1/gamma(1))*y0_real*exp(eps0_combination)/B(1)/lambda1;

    for i=1:N(1)
        x1(i)=max(left*alpha1(i)/(mean_w1(i)/mean_p),0.001);
    end
    
    
end
