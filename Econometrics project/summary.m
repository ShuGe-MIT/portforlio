clear;
%% initial_condition_summary
% 1. parameters
[n, obs, n_cluster, cum_cluster, count_cluster, M, N, A, B,...
    gamma, theta,...
    alpha1, alpha2, alpha3, ratio1, ratio2, ratio3, ratio0]=parameter();
% 2. input and output
load('result_with_regression_little_sai.mat','lnxy_dist','sai_dist_data');
lnxy_dist_sai=[lnxy_dist(:,1:N(1))+squeeze(sai_dist_data(:,1:N(1))),...
                    lnxy_dist(:,(N(1)+1):(N(1)+N(2)))+squeeze(sai_dist_data(:,(N(1)+1):(N(1)+N(2)))),...
                    lnxy_dist(:,(N(1)+N(2)+1):(N(1)+N(2)+N(3)))+squeeze(sai_dist_data(:,(N(1)+N(2)+1):(N(1)+N(2)+N(3)))),...
                    lnxy_dist(:,(N(1)+N(2)+N(3)+1))+squeeze(sai_dist_data(:,(N(1)+N(2)+N(3)+1))),...
                    lnxy_dist(:,(N(1)+N(2)+N(3)+2))+squeeze(sai_dist_data(:,(N(1)+N(2)+N(3)+2))),...
                    lnxy_dist(:,(N(1)+N(2)+N(3)+3))+squeeze(sai_dist_data(:,(N(1)+N(2)+N(3)+3)))];
                
xy_dist=exp(lnxy_dist_sai);
mean(xy_dist)
std(xy_dist)
% 3. Summary for y0 and weather 
load('result_with_regression_little_sai.mat','initial','rain_temp_dist','price_dist');
mean(initial)
std(initial)
mean(rain_temp_dist)
std(rain_temp_dist)
mean(price_dist)
std(price_dist)

%% results
% 1. initial_guess_?
load('guess_valid_new.mat','guess')
data=guess(1:25,1:22);
stats(:,1)=mean(data);
stats(:,2)=std(data);
stats(:,3)=min(data);
stats(:,4)=max(data);
stats(:,5)=quantile(data,0.25);
stats(:,6)=quantile(data,0.5);
stats(:,7)=quantile(data,0.75);

% 2. ?
load('1_5.mat','estimates_beta','error');
% data=estimates_beta((error(1:22)<500),:);
data=estimates_beta(1:5,:);
stats(1:51,1:9)=0.0;
stats(:,1)=mean(data);
stats(:,2)=std(data);
stats(:,3)=min(data);
stats(:,4)=max(data);
stats(:,5)=quantile(data,0.25);
stats(:,6)=quantile(data,0.5);
stats(:,7)=quantile(data,0.75);
stats(:,8)=stats(:,1)./stats(:,2);
stats(:,9)= (1-tcdf(abs(stats(:,8)),MC-1))*2;

% 3. ?
load('1_5.mat','estimates_delta','error');
% data=estimates_delta((error(1:22)<500),:);
data=estimates_delta(1:5,:);
stats(1:204,1:9)=0.0;
stats(:,1)=mean(data);
stats(:,2)=std(data);
stats(:,3)=min(data);
stats(:,4)=max(data);
stats(:,5)=quantile(data,0.25);
stats(:,6)=quantile(data,0.5);
stats(:,7)=quantile(data,0.75);
stats(:,8)=stats(:,1)./stats(:,2);
stats(:,9)= (1-tcdf(abs(stats(:,8)),MC-1))*2;


%% compute stats: mean, std, min, max, 25, 50, 75, t-stat, p-value
stats(1:152,1:9)=0.0;
stats(:,1)=mean(data);
stats(:,2)=std(data);
stats(:,3)=min(data);
stats(:,4)=max(data);
stats(:,5)=quantile(data,0.25);
stats(:,6)=quantile(data,0.5);
stats(:,7)=quantile(data,0.75);
stats(:,8)=stats(:,1)./stats(:,2).*MC^0.5;
stats(:,9)= (1-tcdf(abs(stats(:,8)),MC-1))*2;