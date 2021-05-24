function [x3] = solve_x3(g3,gamma, theta, alpha1, alpha2, alpha3, A, B, ratio, ...
                    eps0_sup, eps0_prob, mean_eps0, ...
                    eps1_sup, eps1_prob, mean_eps1, ...
                    eps2_sup, eps2_prob, mean_eps2, mean_exp_eps3_combination, y2_real, eps0_real, eps1_real, eps2_real, price)
    % this function solves stage 3 input values given structural parameter
    % values
    [n, obs, ~, ~, ~, M, N, ~, ~, ~, ~, ~,...
    ~, ~,...
    ~, ~, ~,~]=parameter();
        
    %mean_w1=price(1:N(1));
    %mean_w2=price((N(1)+1):(N(1)+N(2)));
    mean_w3=price((N(1)+N(2)+1):(N(1)+N(2)+N(3)));
    mean_p=price((N(1)+N(2)+N(3))+1);
    
    x3(1:N(3))=0.0;
    
    lambda3=cal_lambda(N(3), alpha3, price((N(1)+N(2)+1):(N(1)+N(2)+N(3))), price((N(1)+N(2)+N(3))+1));    deviation0=[abs(eps0_real(1)-mean_eps0(1))/mean_eps0(1),...
                abs(eps0_real(2)-mean_eps0(2))/mean_eps0(2)];
    deviation1=[abs(eps1_real(1)-mean_eps1(1))/mean_eps1(1),...
                abs(eps1_real(2)-mean_eps1(2))/mean_eps1(2)];
    deviation2=[abs(eps2_real(1)-mean_eps2(1))/mean_eps2(1),...
                abs(eps2_real(2)-mean_eps2(2))/mean_eps2(2)];
    eps3_combination=1-ratio(2,:)*deviation0'-ratio(3,:)*deviation1'-ratio(4,:)*deviation2';
    eps2_combination=1-ratio(1,:)*deviation2';
    
%     g3=(1-theta(3))^(1/(gamma(3)-1))*(A(3)*B(3)*lambda3*exp(eps3_combination))^(gamma(3)/(gamma(3)-1))-1;
    left=(theta(3)/(1-theta(3))/g3)^(1/gamma(3))*y2_real*exp(eps2_combination)/B(3)/lambda3;
    for i=1:N(3)
        x3(i)=max(left*alpha3(i)/(mean_w3(i)/mean_p),0.001);
    end
end

