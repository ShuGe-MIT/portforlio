% This file is the main file. It calls all the other functions.

%% Housekeeping
clear all;
MC=1;
%% Parameter
% tic;
[n, obs, n_cluster, cum_cluster, count_cluster, M, N, A, B,...
    gamma, theta,...
    alpha1, alpha2, alpha3, ratio1, ratio2, ratio3, ratio0]=parameter();
[w1_sup, w2_sup, w3_sup, p_sup, w1_prob, w2_prob, w3_prob, p_prob]=price(n);
[area_sup, cec_sup, organic_sup, area_prob, cec_prob, organic_prob,...
    eps0_sup, eps1_sup, eps2_sup, eps3_sup, eps0_prob, eps1_prob, eps2_prob, eps3_prob,...
        mean_initial, mean_eps0, mean_eps1, mean_eps2, mean_eps3,...
        mean_exp_eps0_combination, mean_exp_eps1_combination,...
        mean_exp_eps2_combination, mean_exp_eps3_combination]=shock(n);


%% draw initial data
% %% first, draw factor price for each agent xy_dist
% price_dist(1:obs, 1:(sum(N)+1))=0.0;
% price_dist(:,1)=randsrc(1,obs,[squeeze(w1_sup(1,1,:))'; squeeze(w1_prob(1,1,:))'])';
% price_dist(:,2)=randsrc(1,obs,[squeeze(w1_sup(1,2,:))'; squeeze(w1_prob(1,2,:))'])';
% price_dist(:,3)=randsrc(1,obs,[squeeze(w1_sup(1,3,:))'; squeeze(w1_prob(1,3,:))'])';
% price_dist(:,4)=randsrc(1,obs,[squeeze(w1_sup(1,4,:))'; squeeze(w1_prob(1,4,:))'])';
% price_dist(:,5)=randsrc(1,obs,[squeeze(w1_sup(1,5,:))'; squeeze(w1_prob(1,5,:))'])';
% price_dist(:,6)=randsrc(1,obs,[squeeze(w2_sup(1,1,:))'; squeeze(w2_prob(1,1,:))'])';
% price_dist(:,7)=randsrc(1,obs,[squeeze(w2_sup(1,2,:))'; squeeze(w2_prob(1,2,:))'])';
% price_dist(:,8)=randsrc(1,obs,[squeeze(w2_sup(1,3,:))'; squeeze(w2_prob(1,3,:))'])';
% price_dist(:,9)=randsrc(1,obs,[squeeze(w3_sup(1,1,:))'; squeeze(w3_prob(1,1,:))'])';
% price_dist(:,10)=randsrc(1,obs,[squeeze(w3_sup(1,2,:))'; squeeze(w3_prob(1,2,:))'])';
% price_dist(:,11)=randsrc(1,obs,[squeeze(p_sup(1,:)); squeeze(p_prob(1,:))])';
% std(price_dist);
% save('price_data.mat','price_dist');
load('price_data.mat');

% %% second, generate eps for each agent
% rain_temp_index(1:obs, 1:8)=0.0;
% rain_temp_dist(1:obs, 1:8)=0.0;
% rain_temp_index(:,1)=randsrc(1,obs,[1:1:10; squeeze(eps0_prob(1,:))])';
% rain_temp_index(:,2)=randsrc(1,obs,[1:1:10; squeeze(eps0_prob(2,:))])';
% rain_temp_index(:,3)=randsrc(1,obs,[1:1:10; squeeze(eps1_prob(1,:))])';
% rain_temp_index(:,4)=randsrc(1,obs,[1:1:10; squeeze(eps1_prob(2,:))])';
% rain_temp_index(:,5)=randsrc(1,obs,[1:1:10; squeeze(eps2_prob(1,:))])';
% rain_temp_index(:,6)=randsrc(1,obs,[1:1:10; squeeze(eps2_prob(2,:))])';
% rain_temp_index(:,7)=randsrc(1,obs,[1:1:10; squeeze(eps3_prob(1,:))])';
% rain_temp_index(:,8)=randsrc(1,obs,[1:1:10; squeeze(eps3_prob(2,:))])';
% rain_temp_dist(:,1)=eps0_sup(1, rain_temp_index(:,1));
% rain_temp_dist(:,2)=eps0_sup(2, rain_temp_index(:,2));
% rain_temp_dist(:,3)=eps1_sup(1, rain_temp_index(:,3));
% rain_temp_dist(:,4)=eps1_sup(2, rain_temp_index(:,4));
% rain_temp_dist(:,5)=eps2_sup(1, rain_temp_index(:,5));
% rain_temp_dist(:,6)=eps2_sup(2, rain_temp_index(:,6));
% rain_temp_dist(:,7)=eps3_sup(1, rain_temp_index(:,7));
% rain_temp_dist(:,8)=eps3_sup(2, rain_temp_index(:,8));
% std(rain_temp_dist);
% eps_dist(1:obs, 1:8)=0.0;
% for i=1:obs
%     eps_dist(i,1)=abs(rain_temp_dist(i,1)-mean_eps0(1))/mean_eps0(1);
%     eps_dist(i,2)=abs(rain_temp_dist(i,2)-mean_eps0(2))/mean_eps0(2);
%     eps_dist(i,3)=abs(rain_temp_dist(i,3)-mean_eps1(1))/mean_eps1(1);
%     eps_dist(i,4)=abs(rain_temp_dist(i,4)-mean_eps1(2))/mean_eps1(2);
%     eps_dist(i,5)=abs(rain_temp_dist(i,5)-mean_eps2(1))/mean_eps2(1);
%     eps_dist(i,6)=abs(rain_temp_dist(i,6)-mean_eps2(2))/mean_eps2(2);
%     eps_dist(i,7)=abs(rain_temp_dist(i,7)-mean_eps3(1))/mean_eps3(1);
%     eps_dist(i,8)=abs(rain_temp_dist(i,8)-mean_eps3(2))/mean_eps3(2);
% end
% save('eps_data.mat','rain_temp_index','rain_temp_dist','eps_dist');
load('eps_data.mat');

% %% third, generate initial condition
% initial(1:obs,1:3)=0.0;
% initial(:,1)=randsrc(1,obs,[area_sup; area_prob])';
% initial(:,2)=randsrc(1,obs,[cec_sup; cec_prob])';
% initial(:,3)=randsrc(1,obs,[organic_sup; organic_prob])';
% std(initial);
% y0=sum(initial')';
% save('y0.mat','y0','initial');
load('y0.mat');
% 
% x1_dist(1:obs, 1:N(1))=0.0;
% x2_dist(1:obs, 1:N(2))=0.0;
% x3_dist(1:obs, 1:N(3))=0.0;
% y1_dist(1:obs)=0.0;
% y2_dist(1:obs)=0.0;
% y3_dist(1:obs)=0.0;
% 
% toc;
% 
%% solve the model
% tic;
% for i_obs=1:obs
%     y_real0=ratio0*initial(i_obs,:)';
%     rain_temp0=rain_temp_dist(i_obs, 1:2);
%     eps_real0=1-ratio1(1,:)*eps_dist(i_obs,1:2)';
%     
%     % Calculate x1 vector
%     [x1]=solve_x1(gamma, theta, alpha1, alpha2, alpha3, A, B, ratio1, ...
%                     eps0_sup, eps0_prob, mean_eps0, ...
%                     mean_exp_eps1_combination, y_real0, rain_temp0, price_dist(i_obs,:));
%     x1_dist(i_obs,:)=x1;
% 
%     prod_x1=1.0;
%     for i=1:N(1)
%         prod_x1=prod_x1*x1(i)^alpha1(i);
%     end
% 
%     y_real1=A(1)*(theta(1)*(y_real0*exp(eps_real0))^gamma(1)+...
%        (1-theta(1))*(B(1)*prod_x1)^gamma(1))^(1/gamma(1));
%     y1_dist(i_obs)=y_real1;
% 
%     % Calculate x2 vector
%     x2(1:N(2))=0.0;
%     rain_temp1=rain_temp_dist(i_obs, 3:4);
%     eps1_combination=1-ratio2(1,:)*eps_dist(i_obs,3:4)';
%     x2=solve_x2(gamma, theta, alpha1, alpha2, alpha3, A, B, ratio2, ...
%                     eps0_sup, eps0_prob, mean_eps0, ...
%                     eps1_sup, eps1_prob, mean_eps1, ...
%                         mean_exp_eps2_combination, y_real1, rain_temp0, rain_temp1, price_dist(i_obs,:));
% 
%     x2_dist(i_obs,:)=x2;
% 
%     prod_x2=1.0;
%     for i=1:N(2)
%         prod_x2=prod_x2*x2(i)^alpha2(i);
%     end
% 
%     y_real2=0.0;
%     y_real2=A(2)*(theta(2)*(y_real1*exp(eps1_combination))^gamma(2)+...
%             (1-theta(2))*(B(2)*prod_x2)^gamma(2))^(1/gamma(2));
%     y2_dist(i_obs)=y_real2;
% 
%     % Calculate x3 vector
%     x3(1:N(3))=0.0; 
%     y_real3=0.0;
%     rain_temp2=rain_temp_dist(i_obs, 5:6);
%     eps2_combination=1-ratio3(1,:)*eps_dist(i_obs,5:6)';
%     
%     x3=solve_x3(gamma, theta, alpha1, alpha2, alpha3, A, B, ratio3, ...
%                     eps0_sup, eps0_prob, mean_eps0,...
%                     eps1_sup, eps1_prob, mean_eps1,...
%                     eps2_sup, eps2_prob, mean_eps2,...
%                     mean_exp_eps2_combination, y_real2, rain_temp0, rain_temp1, rain_temp2, price_dist(i_obs,:));
%     
% 
%     x3_dist(i_obs,:)=x3;
%     prod_x3=1.0;
%     for i=1:N(3)
%         prod_x3=prod_x3*x3(i)^alpha3(i);
%     end
%     y_real3=A(3)*(theta(3)*(y_real2*exp(eps2_combination))^gamma(3)+...
%             (1-theta(3))*(B(3)*prod_x3)^gamma(3))^(1/gamma(3));
%     y3_dist(i_obs)=y_real3;
% end
% xy_dist=[x1_dist, x2_dist, x3_dist, y1_dist', y2_dist', y3_dist'];
% lnxy_dist=log([x1_dist, x2_dist, x3_dist, y1_dist', y2_dist', y3_dist']);
% %std_lnxy=std(lnxy_dist);
% mu=zeros(sum(N)+3,1);
% %%% construct the cov_mat
% cov_mat=construct_S(n,sum(N)+3);
% index=0.1;
% cov_mat_use=index*cov_mat;
% toc;
%% draw the error terms
% % sai_dist_data(1:obs,1:sum(N)+3)=0.0;
% % sai_dist_data(:,:)=mvnrnd(mu,cov_mat_use,obs);
% % 
% % sai_dist(1:MC,1:M,1:obs,1:sum(N)+3)=0.0;
% % for i=1:MC
% %     for j=1:M
% %         sai_dist(i,j,:,:)=mvnrnd(mu,cov_mat_use,obs);
% %     end
% % end
% %save('result_with_regression_little_sai.mat','-v7.3');
load('result_with_regression_little_sai.mat');

%% preliminaries of matching
no_equations=sum(N)+1;
sai_dist_data=sai_dist_data(:,:);
sai_dist=sai_dist(1:MC,:,:,:);

[n, obs, n_cluster, cum_cluster, count_cluster, M, N, A, B,...
    gamma, theta,...
    alpha1, alpha2, alpha3, ratio1, ratio2, ratio3, ratio0]=parameter();
data=[price_dist, eps_dist, rain_temp_dist, initial, ...
    log(price_dist), log(rain_temp_dist), log(initial)];
save('data.mat','data');

lnxy_dist_sai=[lnxy_dist(:,1:N(1))+squeeze(sai_dist_data(:,1:N(1))),...
                    lnxy_dist(:,(N(1)+1):(N(1)+N(2)))+squeeze(sai_dist_data(:,(N(1)+1):(N(1)+N(2)))),...
                    lnxy_dist(:,(N(1)+N(2)+1):(N(1)+N(2)+N(3)))+squeeze(sai_dist_data(:,(N(1)+N(2)+1):(N(1)+N(2)+N(3)))),...
                    lnxy_dist(:,(N(1)+N(2)+N(3)+1))+squeeze(sai_dist_data(:,(N(1)+N(2)+N(3)+1))),...
                    lnxy_dist(:,(N(1)+N(2)+N(3)+2))+squeeze(sai_dist_data(:,(N(1)+N(2)+N(3)+2))),...
                    lnxy_dist(:,(N(1)+N(2)+N(3)+3))+squeeze(sai_dist_data(:,(N(1)+N(2)+N(3)+3)))];
                
log_price_dist=log(price_dist);
help=[mean_eps0, mean_eps1, mean_eps2, mean_eps3];
log_rain_temp_dist=abs(rain_temp_dist-help)./help;
[BETA_data_pre, W_data_pre]=SUR(no_equations, obs, initial,  log_price_dist, log_rain_temp_dist, lnxy_dist_sai);
initial_coe_guess_x11=BETA_data_pre(120:122)*1/BETA_data_pre(120);
initial_coe_guess_x12=BETA_data_pre(137:139)*1/BETA_data_pre(137);
initial_coe_guess_x13=BETA_data_pre(154:156)*1/BETA_data_pre(154);
initial_coe_guess_x14=BETA_data_pre(171:173)*1/BETA_data_pre(171);
initial_coe_guess_x15=BETA_data_pre(188:190)*1/BETA_data_pre(188);
initial_guess=mean([initial_coe_guess_x11,initial_coe_guess_x12,initial_coe_guess_x13,initial_coe_guess_x14,initial_coe_guess_x15]');

input_choice=mean(exp(lnxy_dist_sai(:,1:sum(N))));
mean_price=mean(price_dist(:,1:sum(N)));
alpha_guess=input_choice.*mean_price;
alpha1_guess=alpha_guess(1:N(1))./sum(alpha_guess(1:N(1)));
alpha2_guess=alpha_guess((N(1)+1):(N(1)+N(2)))./sum(alpha_guess((N(1)+1):(N(1)+N(2))));
alpha3_guess=alpha_guess((N(1)+N(2)+1):(N(1)+N(2)+N(3)))./sum(alpha_guess((N(1)+N(2)+1):(N(1)+N(2)+N(3))));
% save('guess_initial_alpha.mat','initial_guess','alpha1_guess','alpha2_guess','alpha3_guess');

log_initial=log(initial);
log_price_dist=log(price_dist);
log_rain_temp_dist=log(rain_temp_dist);
% calculate beta and covariance matrix with true structural parameter values
[BETA_data, COV_data]=SUR_real_matrix(no_equations, obs, log_initial,  log_price_dist, log_rain_temp_dist, lnxy_dist_sai);
W_data=inv(COV_data);
% W_data*COV_data;
%[BETA_data, COV_data]=SUR_synthetic(no_equations, obs, log_initial,  log_price_dist, log_rain_temp_dist, lnxy_dist_sai);
%[BETA_data, W_data]=SUR(no_equations, obs, log_initial,  log_price_dist, log_rain_temp_dist, lnxy_dist_sai);

estimates_beta(1:MC,1:51)=0.0;
estimates_delta(1:MC,1:204)=0.0;
error(1:MC)=0.0;
exitflag(1:MC)=0.0;
funeval(1:MC)=0.0;
% para_guess=[gamma, theta, alpha1, alpha2, alpha3, A, B,...
%     ratio1(1,:), ratio1(2,:), ratio1(3,:), ratio1(4,:),...
%     ratio2(1,:), ratio2(2,:), ratio2(3,:), ratio2(4,:), ratio2(5,:),...
%     ratio3(1,:), ratio3(2,:), ratio3(3,:), ratio3(4,:), ratio0];
% compute_error(para_guess, temp, data, BETA_data, W_data);
test_start=1;
test_end=1;
initial_error(1:test_end)=0.0;
load('guess_polished.mat','guess_polished');

%% matching
% This section tries to find the structural parameter values
for i=test_start:test_end
    i    
    tic;
    lb=[-10,-10,-10,0,0,0,zeros(1,sum(N)),0.5,0.5,0.5,0.5,0.5,0.5,ones(1,26)*(-3),0,0,0];
    ub=[1,1,1,1,1,1,ones(1,sum(N)),5.0,5.0,5.0,5.0,5.0,5.0,ones(1,26)*(3),5,5,5];
    temp=squeeze(sai_dist(i,:,:,:));
    fun = @(para) compute_error(para, temp, data, BETA_data, W_data);
    options = optimset('MaxFunEvals',5000,'Algorithm','interior-point','OutputFcn',@outfun);
    para_guess=guess_polished(i,:);
    initial_error(i)=compute_error(para_guess, temp, data, BETA_data, W_data);
	[estimates_beta(i,:),error(i),exitflag(i),c]= fmincon(fun, para_guess, [],[],[],[],lb,ub,[],options);
    funeval(i)=c.funcCount;
    toc;

end

% After finding the structural parameter values, calculate the beta (auxiliary parameter values)
BETA_cal(1:M,1:204)=0.0;
for i=test_start:test_end
    [xy_dist_cal]=solve_x(estimates_beta(i,:), data(:,1:(sum(N)+20)));
    lnxy_dist_cal=log(xy_dist_cal);
    lnxy_dist_cal_sai(1:M,1:obs,1:sum(N)+3)=0.0;
    for j=1:M
        lnxy_dist_cal_sai(j,:,:)=lnxy_dist_cal+squeeze(sai_dist(i,j,:,:));
        temp=squeeze(lnxy_dist_cal_sai(j,:,:));
        [BETA_cal(j,:), ~]=SUR(no_equations, obs, log_initial, log_price_dist, log_rain_temp_dist, temp);
        % on the above line, it does not matter whether it is SUR, or
        % SUR_real_matrix, or SUR_synthetic since we only want BETA
    end
    estimates_delta(i,:)=mean(BETA_cal);
end
toc;

save('1_5.mat','estimates_beta','estimates_delta','error','exitflag','funeval');
save('1_5_complete.mat');
