function [xy_dist_cal] = solve_x(para, data)
    % this function solves the model with given structural parameter values
    % and given realized shocks.
    [n, obs, ~, ~, ~, M, N, ~, ~, ~, ~, ~,...
    ~, ~,...
    ~, ~, ~, ~]=parameter();

    [area_sup, cec_sup, organic_sup, area_prob, cec_prob, organic_prob,...
     eps0_sup, eps1_sup, eps2_sup, eps3_sup, eps0_prob, eps1_prob, eps2_prob, eps3_prob,...
        mean_initial, mean_eps0, mean_eps1, mean_eps2, mean_eps3,...
        mean_exp_eps0_combination, mean_exp_eps1_combination,...
        mean_exp_eps2_combination, mean_exp_eps3_combination]=shock(n);
    
    price_dist=data(:,1:(sum(N)+1));
    eps_dist=data(:,(sum(N)+2):(sum(N)+9));
    rain_temp_dist=data(:,(sum(N)+10):(sum(N)+17));
    initial=data(:,(sum(N)+18):(sum(N)+20));

    x1_dist_cal(1:obs, 1:N(1))=0.0;
    x2_dist_cal(1:obs, 1:N(2))=0.0;
    x3_dist_cal(1:obs, 1:N(3))=0.0;
    y1_dist_cal(1:obs)=0.0;
    y2_dist_cal(1:obs)=0.0;
    y3_dist_cal(1:obs)=0.0;

    gamma=para(1:3);
    theta=para(4:6);
    alpha1=para(7:(N(1)+6));
    alpha2=para((N(1)+7):(N(1)+N(2)+6));
    alpha3=para((N(1)+N(2)+7):(N(1)+N(2)+N(3)+6));
    A=para((N(1)+N(2)+N(3)+7):(N(1)+N(2)+N(3)+9));
    B=para((N(1)+N(2)+N(3)+10):(N(1)+N(2)+N(3)+12));
    help=sum(N)+12;
    ratio1=[para(help+[1,3,5,7])',para(help+[2,4,6,8])'];
    ratio2=[para(help+[9,11,13,15,17])',para(help+[10,12,14,16,18])'];
    ratio3=[para(help+[19,21,23,25])',para(help+[20,22,24,26])'];
    ratio0=para(help+26+[1,2,3]);
        
    for i_obs=1:obs
        %% calculate x1
        y_real0=ratio0*initial(i_obs,:)';
        rain_temp0=rain_temp_dist(i_obs, 1:2);
        eps_real0=1-ratio1(1,:)*eps_dist(i_obs,1:2)';

        % Calculate x1 vector
        [x1_cal]=solve_x1_old(gamma, theta, alpha1, alpha2, alpha3, A, B, ratio1, ...
                        eps0_sup, eps0_prob, mean_eps0, ...
                    mean_exp_eps1_combination, y_real0, rain_temp0, price_dist(i_obs,:));
        
%         if (x1_cal~=real(x1_cal))
%             x1_cal(1:N(1))=0.0;
%         end
        x1_dist_cal(i_obs,:)=x1_cal;
        
        prod_x1=1.0;
        for i=1:N(1)
            prod_x1=prod_x1*x1_cal(i)^alpha1(i);
        end

        y_real1=A(1)*(theta(1)*(y_real0*exp(eps_real0))^gamma(1)+...
           (1-theta(1))*(B(1)*prod_x1)^gamma(1))^(1/gamma(1));
        y1_dist_cal(i_obs)=y_real1;
        %% calculate x2
        x2_cal(1:N(2))=0.0;
        rain_temp1=rain_temp_dist(i_obs, 3:4);
        eps1_combination=1-ratio2(1,:)*eps_dist(i_obs,3:4)';
        x2_cal=solve_x2_old(gamma, theta, alpha1, alpha2, alpha3, A, B, ratio2, ...
                        eps0_sup, eps0_prob, mean_eps0, ...
                        eps1_sup, eps1_prob, mean_eps1, ...
                        mean_exp_eps2_combination, y_real1, rain_temp0, rain_temp1, price_dist(i_obs,:));
        
%         if (x2_cal~=real(x2_cal))
%             x2_cal(1:N(2))=0.0;
%         end
        x2_dist_cal(i_obs,:)=x2_cal;

        prod_x2=1.0;
        for i=1:N(2)
            prod_x2=prod_x2*x2_cal(i)^alpha2(i);
        end

        y_real2=0.0;
        y_real2=A(2)*(theta(2)*(y_real1*exp(eps1_combination))^gamma(2)+...
                (1-theta(2))*(B(2)*prod_x2)^gamma(2))^(1/gamma(2));

        y2_dist_cal(i_obs)=y_real2;
        %% Calculate x3 vector
        x3_cal(1:N(3))=0.0; 
        y_real3=0.0;
        rain_temp2=rain_temp_dist(i_obs, 5:6);
        eps2_combination=1-ratio3(1,:)*eps_dist(i_obs,5:6)';
        x3_cal=solve_x3_old(gamma, theta, alpha1, alpha2, alpha3, A, B, ratio3, ...
                    eps0_sup, eps0_prob, mean_eps0,...
                    eps1_sup, eps1_prob, mean_eps1,...
                    eps2_sup, eps2_prob, mean_eps2,...
                mean_exp_eps2_combination, y_real2, rain_temp0, rain_temp1, rain_temp2, price_dist(i_obs,:));
        prod_x3=1.0;
        for i=1:N(3)
            prod_x3=prod_x3*x3_cal(i)^alpha3(i);
        end
%         if (x3_cal~=real(x3_cal))
%             x3_cal(1:N(3))=0.0;
%         end
        y_real3=A(3)*(theta(3)*(y_real2*exp(eps2_combination))^gamma(3)+...
                (1-theta(3))*(B(3)*prod_x3)^gamma(3))^(1/gamma(3));

        x3_dist_cal(i_obs,:)=x3_cal;
        y3_dist_cal(i_obs)=y_real3;
    end
%     alpha1;
%     alpha2;
%     alpha3;
%     theta;
%     gamma;
%     A;
%     B;
%     ratio1;
%     ratio2;
%     ratio3;
    xy_dist_cal=[x1_dist_cal, x2_dist_cal, x3_dist_cal, y1_dist_cal', y2_dist_cal', y3_dist_cal'];

end

