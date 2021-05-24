function [Wald,G] = compute_Wald(para)
    % ,sai_dist,data,DELTA_data,V_data
    % para is the guess of structural parameter values
    % sai_dist is the Wald term
    % data are all the random terms draVn in the beginning, including prices, shocks, initial
    % conditions
    load('V_data.mat','V_data');
    load('tracker','iteration','INITIAL_DIFF','INITIAL','INITIAL_UPDATE','DELTA_BAR','DELTA_DIFF','W_sumsq','W_eig','W_sumabs','WALD');
    load('DELTA_data','DELTA_data');
    load('data.mat','data');
    load('sai_dist.mat','temp');
    sai_dist=temp;
    load('true_para');
    save('current_para','para');
    iteration=iteration+1;
    INITIAL(:,iteration)=para;
%%    
    [n, obs, ~, ~, ~, M, N, ~, ~, ~, ~, ~,...
    ~, ~,...
    ~, ~, ~, ~]=parameter();
    log_price_dist=data(:,(sum(N)+21):(2*sum(N)+21));
    log_rain_temp=data(:,(2*sum(N)+22):(2*sum(N)+29));
    log_initial=data(:,(2*sum(N)+30):(2*sum(N)+32));
    % calculate input and output values
    data=data(:,1:(sum(N)+20));
    % this function solves the model with given structural parameter values
    % and given realized shocks.
    [area_sup, cec_sup, organic_sup, area_prob, cec_prob, organic_prob,...
     eps0_sup, eps1_sup, eps2_sup, eps3_sup, eps0_prob, eps1_prob, eps2_prob, eps3_prob,...
        mean_initial, mean_eps0, mean_eps1, mean_eps2, mean_eps3,...
        mean_exp_eps0_combination, mean_exp_eps1_combination,...
        mean_exp_eps2_combination, mean_exp_eps3_combination]=shock(n);
    price_dist=data(:,1:(sum(N)+1));
    eps_dist=data(:,(sum(N)+2):(sum(N)+9));
    rain_temp_dist=data(:,(sum(N)+10):(sum(N)+17));
    initial=data(:,(sum(N)+18):(sum(N)+20));
% save('test_dist_cal','price_dist','eps_dist','rain_temp_dist','initial');
    x1_dist_cal(1:obs, 1:N(1))=0.0;
    x2_dist_cal(1:obs, 1:N(2))=0.0;
    x3_dist_cal(1:obs, 1:N(3))=0.0;
    y1_dist_cal(1:obs)=0.0;
    y2_dist_cal(1:obs)=0.0;
    y3_dist_cal(1:obs)=0.0;
    gamma=para(1:3)';
    theta=para(4:6)';
    alpha1=para(7:(N(1)+6))';
    alpha2=para((N(1)+7):(N(1)+N(2)+6))';
    alpha3=para((N(1)+N(2)+7):(N(1)+N(2)+N(3)+6))';
    A=para((N(1)+N(2)+N(3)+7):(N(1)+N(2)+N(3)+9))';
    B=para((N(1)+N(2)+N(3)+10):(N(1)+N(2)+N(3)+12))';
    help=sum(N)+12;
    ratio1=[para(help+[1,3,5,7]),para(help+[2,4,6,8])];
    ratio2=[para(help+[9,11,13,15,17]),para(help+[10,12,14,16,18])];
    ratio3=[para(help+[19,21,23,25]),para(help+[20,22,24,26])];
    ratio0=para(help+26+[1,2,3])';
%% theoretical restrictions
    margin=0.02;
    gamma=min(gamma,1-margin);%gamma <=1
    gamma(abs(gamma)<margin)=margin;
    theta=min(theta,1-margin);%theta <1
    theta=max(theta,margin);%theta >0
    alpha1=min(alpha1,1-margin);%alpha <1
    alpha1=max(alpha1,margin);%alpha >0
    alpha2=min(alpha2,1-margin);%alpha <1
    alpha2=max(alpha2,margin);%alpha >0
    alpha3=min(alpha3,1-margin);%alpha <1
    alpha3=max(alpha3,margin);%alpha >0
    A=max(A,margin);
    B=max(B,margin);
%% check conditions
target=0.95; %target<1
deviation0=(abs(rain_temp_dist(:,1:2)-ones(2000,1)*mean_eps0(1,1:2)))*[1/mean_eps0(1),0;0,1/mean_eps0(2)];
deviation1=(abs(rain_temp_dist(:,3:4)-ones(2000,1)*mean_eps1(1,1:2)))*[1/mean_eps1(1),0;0,1/mean_eps1(2)];
deviation2=(abs(rain_temp_dist(:,5:6)-ones(2000,1)*mean_eps2(1,1:2)))*[1/mean_eps2(1),0;0,1/mean_eps2(2)];
lambda1=ones(2000,1);
lambda2=ones(2000,1);
lambda3=ones(2000,1);
for i=1:2000
lambda1(i)=cal_lambda(N(1),alpha1,price_dist(i,1:N(1)),price_dist(i,(N(1)+N(2)+N(3)+1)));
lambda2(i)=cal_lambda(N(2), alpha2, price_dist(i,((N(1)+1):(N(1)+N(2)))), price_dist(i,(N(1)+N(2)+N(3))+1));
lambda3(i)=cal_lambda(N(3), alpha3, price_dist(i,(N(1)+N(2)+1):(N(1)+N(2)+N(3))), price_dist(i,(N(1)+N(2)+N(3))+1));
end
eps3_combination3=ones(2000,1)-(ratio3(2,:)*deviation0')'-(ratio3(3,:)*deviation1')'-(ratio3(4,:)*deviation2')';
eps3_combination2=ones(2000,1)-(ratio2(4,:)*deviation0')'-(ratio2(5,:)*deviation1')';
eps2_combination2=ones(2000,1)-(ratio2(2,:)*deviation0')'-(ratio2(3,:)*deviation1')';
eps3_combination1=ones(2000,1)-(ratio1(4,:)*deviation0')';
eps2_combination1=ones(2000,1)-(ratio1(3,:)*deviation0')';
eps1_combination1=ones(2000,1)-(ratio1(2,:)*deviation0')';
%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base3=A(3)*B(3)*lambda3.*exp(eps3_combination3);
G3=(1-theta(3))*((base3).^gamma(3));
[x,y]=max(G3);
if x>=1
% % %     fprintf("gamma3/theta3 is changed\n")
    theta(3)=min(1+(1/log(base3(y))),1-margin);
    gamma(3)=min((log(target)-log(1-theta(3)))/(log(base3(y))),1-margin);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base3=A(3)*B(3)*lambda3.*exp(eps3_combination2);
G3_2=(1-theta(3))*((base3).^gamma(3));
[x,y]=max(G3_2);
if x>=1
% % %     fprintf("gamma3/theta3 is changed\n")
    theta(3)=min(1+(1/log(base3(y))),1-margin);
    gamma(3)=min((log(target)-log(1-theta(3)))/(log(base3(y))),1-margin);
end
value=((A(3)*exp(eps3_combination2)).^(gamma(3)/(gamma(3)-1))-...
        ((1-theta(3))*(B(3)*lambda3).^gamma(3)).^(1/(1-gamma(3)))).^((gamma(3)-1)/gamma(3));   
P3=theta(3)^(1/gamma(3))*value;
base2=A(2)*B(2)*lambda2.*P3.*exp(eps2_combination2);
G2=(1-theta(2))*(base2).^gamma(2);
[x,y]=max(G2);
if x>=1
% % %     fprintf("gamma2/theta2 is changed\n")
    theta(2)=min(1+(1/log(base2(y))),1-margin);
    gamma(2)=min((log(target)-log(1-theta(2)))/(log(base2(y))),1-margin);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base3=A(3)*B(3)*lambda3.*exp(eps3_combination1);
G3_1=(1-theta(3))*((base3).^gamma(3));
[x,y]=max(G3_1);
if x>=1
% % %     fprintf("gamma3/theta3 is changed\n")
    theta(3)=min(1+(1/log(base3(y))),1-margin);
    gamma(3)=min((log(target)-log(1-theta(3)))/(log(base3(y))),1-margin);
end
value1=((A(3)*exp(eps3_combination1)).^(gamma(3)/(gamma(3)-1))-...
           ((1-theta(3))*(B(3)*lambda3).^gamma(3)).^(1/(1-gamma(3)))).^((gamma(3)-1)/gamma(3));
P3_1=theta(3)^(1/gamma(3))*value1;
base2=A(2)*B(2)*lambda2.*P3_1.*exp(eps2_combination1);
G2_1=(1-theta(2))*(base2).^gamma(2);
[x,y]=max(G2_1);
if x>=1
% % %     fprintf("gamma2/theta2 is changed\n")
    theta(2)=min(1+(1/log(base2(y))),1-margin);
    gamma(2)=min((log(target)-log(1-theta(2)))/(log(base2(y))),1-margin);
end
value2=((A(2)*exp(eps2_combination1).*P3_1).^(gamma(2)/(gamma(2)-1))-...
           ((1-theta(2))*(B(2)*lambda2).^gamma(2)).^(1/(1-gamma(2)))).^((gamma(2)-1)/gamma(2));
P2=theta(2)^(1/gamma(2))*value2;
base1=A(1)*B(1)*lambda1.*P2.*exp(eps1_combination1);
G1=(1-theta(1))*(base1).^gamma(1);
[x,y]=max(G1);
if x>=1
% % %     fprintf("gamma1/theta1 is changed\n")
    theta(1)=min(1+(1/log(base1(y))),1-margin);
    gamma(1)=min((log(target)-log(1-theta(1)))/(log(base1(y))),1-margin);
end
%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base3=A(3)*B(3)*lambda3.*exp(eps3_combination3);
G3=(1-theta(3))*(base3).^gamma(3);

value=((A(3)*exp(eps3_combination2)).^(gamma(3)/(gamma(3)-1))-...
        ((1-theta(3))*(B(3)*lambda3).^gamma(3)).^(1/(1-gamma(3)))).^((gamma(3)-1)/gamma(3));    
P3=theta(3)^(1/gamma(3))*value;
base2=A(2)*B(2)*lambda2.*P3.*exp(eps2_combination2);
G2=(1-theta(2))*(base2).^gamma(2);
value1=((A(3)*exp(eps3_combination1)).^(gamma(3)/(gamma(3)-1))-...
           ((1-theta(3))*(B(3)*lambda3).^gamma(3)).^(1/(1-gamma(3)))).^((gamma(3)-1)/gamma(3));
P3_1=theta(3)^(1/gamma(3))*value1;
value2=((A(2)*exp(eps2_combination1).*P3_1).^(gamma(2)/(gamma(2)-1))-...
           ((1-theta(2))*(B(2)*lambda2).^gamma(2)).^(1/(1-gamma(2)))).^((gamma(2)-1)/gamma(2));       
P2=theta(2)^(1/gamma(2))*value2;
base1=A(1)*B(1)*lambda1.*P2.*exp(eps1_combination1);
G1=(1-theta(1))*(base1).^gamma(1);
% % % if G1>1   
% % %     fprintf("there is bug in G1\n")
% % % elseif G2>1
% % %     fprintf("there is bug in G2\n")
% % % elseif G3>1
% % %     fprintf("there is bug in G3\n")
% % % else
% % %     fprintf("there is no bug\n")
% % % end

INITIAL_UPDATE(:,iteration)=[gamma,theta,alpha1,alpha2,alpha3,A,B,para(23:51)']';
INITIAL_DIFF(:,iteration)=abs(true_para-[gamma,theta,alpha1,alpha2,alpha3,A,B,para(23:51)']');
    %% continue calculation        
    for i_obs=1:obs
        %% calculate x1
        y_real0=ratio0*initial(i_obs,:)';
        if y_real0<0
            [~,y]=min(ratio0);
            y_real0(y)=1;
            y_real0=ratio0*initial(i_obs,:)';
        end
        rain_temp0=rain_temp_dist(i_obs, 1:2);
        eps_real0=1-ratio1(1,:)*eps_dist(i_obs,1:2)';
        % Calculate x1 vector
        [x1_cal]=solve_x1(G1(i_obs)^(1/(gamma(1)-1))-1,gamma, theta, alpha1, alpha2, alpha3, A, B, ratio1, ...
                        eps0_sup, eps0_prob, mean_eps0, ...
                    mean_exp_eps1_combination, y_real0, rain_temp0, price_dist(i_obs,:));
%        [x1_cal]=solve_x1_old(gamma, theta, alpha1, alpha2, alpha3, A, B, ratio1, ...
%                     eps0_sup, eps0_prob, mean_eps0, ...
%                     mean_exp_eps1_combination, y_real0, rain_temp0, price_dist(i_obs,:));
% save('test_cal','x1_cal','gamma', 'theta', 'alpha1', 'alpha2', 'alpha3', 'A', 'B', 'ratio1', ...
%                     'eps0_sup', 'eps0_prob', 'mean_eps0', ...
%                     'mean_exp_eps1_combination', 'y_real0', 'rain_temp0', 'price_dist');
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
        x2_cal=solve_x2(G2(i_obs)^(1/(gamma(2)-1))-1,gamma, theta, alpha1, alpha2, alpha3, A, B, ratio2, ...
                        eps0_sup, eps0_prob, mean_eps0, ...
                        eps1_sup, eps1_prob, mean_eps1, ...
                        mean_exp_eps2_combination, y_real1, rain_temp0, rain_temp1, price_dist(i_obs,:));
        
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
%         g3=G3.^(1/(gamma(3)-1))-1;
%         save('g3','g3');
        x3_cal(1:N(3))=0.0; 
        y_real3=0.0;
        rain_temp2=rain_temp_dist(i_obs, 5:6);
        eps2_combination=1-ratio3(1,:)*eps_dist(i_obs,5:6)';
        x3_cal=solve_x3(G3(i_obs)^(1/(gamma(3)-1))-1,gamma, theta, alpha1, alpha2, alpha3, A, B, ratio3, ...
                    eps0_sup, eps0_prob, mean_eps0,...
                    eps1_sup, eps1_prob, mean_eps1,...
                    eps2_sup, eps2_prob, mean_eps2,...
                mean_exp_eps2_combination, y_real2, rain_temp0, rain_temp1, rain_temp2, price_dist(i_obs,:));
        prod_x3=1.0;
        y_real3=A(3)*(theta(3)*(y_real2*exp(eps2_combination))^gamma(3)+...
                (1-theta(3))*(B(3)*prod_x3)^gamma(3))^(1/gamma(3));

        x3_dist_cal(i_obs,:)=x3_cal;
        y3_dist_cal(i_obs)=y_real3;    
    end
%     save('x3_dist_cal','x3_dist_cal');
    xy_dist_cal=[x1_dist_cal, x2_dist_cal, x3_dist_cal, y1_dist_cal', y2_dist_cal', y3_dist_cal'];
    save('results_synthetic','xy_dist_cal');
    save('results_x3_dist_cal','x3_dist_cal');
%     [xy_dist_cal]=solve_x(para, data(:,1:(sum(N)+20)));
    lnxy_dist_cal=log(xy_dist_cal);
    lnxy_dist_cal_sai(1:M,1:obs,1:(sum(N)+3))=0.0;
    DELTA_cal(1:M,1:204)=0.0;
    V_cal(1:M,1:204,1:204)=0.0;
    no_equations=sum(N)+1;
    parfor i=1:M
        % add Wald terms
        lnxy_dist_cal_sai(i,:,:)=lnxy_dist_cal+squeeze(sai_dist(i,:,:));
        temp=squeeze(lnxy_dist_cal_sai(i,:,:));
        % get DELTA and covariance matrix of each synthetic data
        [DELTA_cal(i,:), COV_cal]=SUR_synthetic(no_equations, obs, log_initial, log_price_dist, log_rain_temp, temp);
        V_cal(i,:,:)=inv(COV_cal);
    end
    V_sum=squeeze(sum(V_cal))/M^2;
    V=V_data+V_sum;
    DELTA_cal_bar=mean(DELTA_cal);
    DELTA_delta=(DELTA_data-DELTA_cal_bar');
    W=inv(V);
    Wald=DELTA_delta'*W*DELTA_delta;
    G=[];
    %save data    
    DELTA_cal_bar_trans=DELTA_cal_bar';
    DELTA_BAR(:,iteration)=DELTA_cal_bar_trans;
    DELTA_DIFF(:,iteration)=DELTA_delta;
    WALD(:,iteration)=real(Wald);
    W_sumsq(:,iteration)=sum(W.^2)';
    W_eig(:,iteration)=eig(W);
    W_sumabs(:,iteration)=sum(abs(W));
    save('tracker','iteration','INITIAL_DIFF','INITIAL','INITIAL_UPDATE','DELTA_BAR','DELTA_DIFF','W_sumsq','W_eig','W_sumabs','WALD');
   end
% %%
% i_obs=0;
%     while i_obs<obs
%         i_obs=i_obs+1;
%         rain_temp0=rain_temp_dist(i_obs, 1:2);
%         rain_temp1=rain_temp_dist(i_obs, 3:4);
%         rain_temp2=rain_temp_dist(i_obs, 5:6);
%         price=price_dist(i_obs,:);
% %% check conditions for g3>0, if not, set base3=0.95
%         lambda1=cal_lambda(N(1), alpha1, price(1:N(1)), price((N(1)+N(2)+N(3))+1));
%         lambda2=cal_lambda(N(2), alpha2, price((N(1)+1):(N(1)+N(2))), price((N(1)+N(2)+N(3))+1));
%         lambda3=cal_lambda(N(3), alpha3, price((N(1)+N(2)+1):(N(1)+N(2)+N(3))), price((N(1)+N(2)+N(3))+1));     
%         deviation0=[abs(rain_temp0(1)-mean_eps0(1))/mean_eps0(1),...
%                 abs(rain_temp0(2)-mean_eps0(2))/mean_eps0(2)];
%         deviation1=[abs(rain_temp1(1)-mean_eps1(1))/mean_eps1(1),...
%                 abs(rain_temp1(2)-mean_eps1(2))/mean_eps1(2)];
%         deviation2=[abs(rain_temp2(1)-mean_eps2(1))/mean_eps2(1),...
%                 abs(rain_temp2(2)-mean_eps2(2))/mean_eps2(2)];
%         eps3_combination=1-ratio3(2,:)*deviation0'-ratio3(3,:)*deviation1'-ratio3(4,:)*deviation2';
%         base3=A(3)*B(3)*lambda3*exp(eps3_combination);
%     if (1-theta(3))*((base3)^gamma(3))>1
%         theta(3)=1+(1/log(base3));
%         gamma(3)=(log(0.95)-log(-1/log(base3)))/(log(base3));
%         i_obs=0; %check from the first observation again
%     end       
% %% check conditions for g2>0, if not, set base2=0.95
%     EP3=0.0; %% Don't forget to time theta3^(1/gamma3)
%     value=0.0; %% store intermediate value for Q3^((gamma3-1)/gamma3)
%     
%     eps3_combination2=1-ratio2(4,:)*deviation0'-ratio2(5,:)*deviation1';
%     eps2_combination2=1-ratio2(2,:)*deviation0'-ratio2(3,:)*deviation1';
%     value=((A(3)*exp(eps3_combination2))^(gamma(3)/(gamma(3)-1))-...
%         ((1-theta(3))*(B(3)*lambda3)^gamma(3))^(1/(1-gamma(3))))^((gamma(3)-1)/gamma(3));
%     EP3=value*theta(3)^(1/gamma(3));
%     base2=A(2)*B(2)*lambda2*EP3*exp(eps2_combination2);
%     if (1-theta(2))*(base2)^gamma(2)>1
%         theta(2)=1+(1/log(base2));
%         gamma(2)=(log(0.95)-log(-1/log(base2)))/(log(base2));
%         i_obs=0;
%     end    
% %% check conditions for g1>0, if not, set base1=0.95
%     EP2=0.0; %% Don't forget to time theta2^(1/gamma2)
%     value3=0.0; %% store intermediate value for Q3^((gamma3-1)/gamma3)
%     value2=0.0; %% stode intermediate value for Q2^((gamma2-1)/gamma2)
%     eps3_combination1=1-ratio1(4,:)*deviation0';
%     eps2_combination1=1-ratio1(3,:)*deviation0';
%     eps1_combination1=1-ratio1(2,:)*deviation0';
%     value3=((A(3)*exp(eps3_combination1))^(gamma(3)/(gamma(3)-1))-...
%            ((1-theta(3))*(B(3)*lambda3)^gamma(3))^(1/(1-gamma(3))))^((gamma(3)-1)/gamma(3));
%     P3=theta(3)^(1/gamma(3))*value3;
%     value2=((A(2)*exp(eps2_combination1)*P3)^(gamma(2)/(gamma(2)-1))-...
%            ((1-theta(2))*(B(2)*lambda2)^gamma(2))^(1/(1-gamma(2))))^((gamma(2)-1)/gamma(2));
%     EP2=value2*theta(2)^(1/gamma(2));
%     base1=A(1)*B(1)*lambda1*EP2*exp(eps1_combination1);
%     if (1-theta(1))*(base1)^gamma(1)>1
%         theta(1)=1+(1/log(base1));
%         gamma(1)=(log(0.95)-log(-1/log(base1)))/(log(base1));  
%         i_obs=0;
%     end   
% end    
