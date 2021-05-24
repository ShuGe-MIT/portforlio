function [x2] = solve_x2(g2,gamma, theta, alpha1, alpha2, alpha3, A, B, ratio, ...
                        eps0_sup, eps0_prob, mean_eps0, ...
                        eps1_sup, eps1_prob, mean_eps1, ...
                        mean_exp_eps2_combination, y1_real, eps0_real, eps1_real, price)
    % this function solves stage 2 input values given structural parameter
    % values
    [n, obs, ~, ~, ~, M, N, ~, ~, ~, ~, ~,...
    ~, ~,...
    ~, ~, ~,~]=parameter();

    %mean_w1=price(1:N(1));
    mean_w2=price((N(1)+1):(N(1)+N(2)));
    mean_w3=price((N(1)+N(2)+1):(N(1)+N(2)+N(3)));
    mean_p=price((N(1)+N(2)+N(3))+1);
    
    x2(1:N(2))=0.0;
    
    lambda2=cal_lambda(N(2), alpha2, price((N(1)+1):(N(1)+N(2))), price((N(1)+N(2)+N(3))+1));
    lambda3=cal_lambda(N(3), alpha3, price((N(1)+N(2)+1):(N(1)+N(2)+N(3))), price((N(1)+N(2)+N(3))+1));
    
    deviation0=[abs(eps0_real(1)-mean_eps0(1))/mean_eps0(1),...
                abs(eps0_real(2)-mean_eps0(2))/mean_eps0(2)];
    deviation1=[abs(eps1_real(1)-mean_eps1(1))/mean_eps1(1),...
                abs(eps1_real(2)-mean_eps1(2))/mean_eps1(2)];
    
    EP3=0.0; %% Don't forget to time theta3^(1/gamma3)
    value=0.0; %% store intermediate value for Q3^((gamma3-1)/gamma3)
    
    eps3_combination=1-ratio(4,:)*deviation0'-ratio(5,:)*deviation1';
    eps2_combination=1-ratio(2,:)*deviation0'-ratio(3,:)*deviation1';
    eps1_combination=1-ratio(1,:)*deviation1';

%     value=((A(3)*exp(eps3_combination))^(gamma(3)/(gamma(3)-1))-...
%         ((1-theta(3))*(B(3)*lambda3)^gamma(3))^(1/(1-gamma(3))))^((gamma(3)-1)/gamma(3));
% 
%     EP3=value*theta(3)^(1/gamma(3));
%     g2=(1-theta(2))^(1/(gamma(2)-1))*...
%         (A(2)*B(2)*lambda2*EP3*exp(eps2_combination))^(gamma(2)/(gamma(2)-1))-1;        
    left=(theta(2)/(1-theta(2))/g2)^(1/gamma(2))*y1_real*exp(eps1_combination)/B(2)/lambda2;
    for i=1:N(2)
        x2(i)=max(left*alpha2(i)/(mean_w2(i)/mean_p),0.001); 
    end

end

