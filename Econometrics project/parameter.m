function [n, obs, n_cluster, cum_cluster, count_cluster, M, N, A, B,...
    gamma, theta,...
    alpha1, alpha2, alpha3, ratio1, ratio2, ratio3, ratio0] = parameter()
    
    % This part draws random structural parameter values, I only draw them
    % once and then save them. 
%     gamma(1:3)=random_num(1,3,0.1,1);
%     
%     theta(1:3)=random_num(1,3,0.2,0.8);
%     
%     alpha3(1:2)=0.0;
%     alpha3(1)=random_num(1,1,0.2,0.7);
%     alpha3(2)=1-alpha3(1);
%     
%     alpha2(1:3)=0.0;
%     alpha2(1)=random_num(1,1,0.2,0.7);
%     alpha2(2)=random_num(1,1,0.0,1-alpha2(1));
%     alpha2(3)=1-alpha2(1)-alpha2(2);
%     
%     alpha1(1:4)=0.0;
%     alpha1(1)=random_num(1,1,0.2,0.7);
%     alpha1(2)=random_num(1,1,0.0,1-alpha1(1));
%     alpha1(3)=random_num(1,1,0.0,1-alpha1(1)-alpha1(2));
%     alpha1(4)=random_num(1,1,0.0,1-alpha1(1)-alpha1(2)-alpha1(3));
%     alpha1(5)=1-alpha1(1)-alpha1(2)-alpha1(3)-alpha1(4);
%     
%     save('structural_parameter.mat');
    
    n=10; % number of discretized grids of each distribution
    obs=2000; % number of observations
    obs_in_each_cluster=4; % number of observations in each cluster
    n_cluster=obs_in_each_cluster.*ones(1,obs/obs_in_each_cluster); % number of clusters
    cum_cluster=cumsum(n_cluster); % total number of clusters
    count_cluster=[250, 250]; % number of clusters in each strata
%     n_cluster=[598, 534, 434, 434];
%     cum_cluster=cumsum(n_cluster);
%     count_cluster=[2,2];
%     n_cluster=[100,300,100,300,100,300,100,300,100,300];
%     cum_cluster=cumsum(n_cluster);
%     count_cluster=[2,3,3,2];
    M=25; % number of sets of synthetic data
    
    N=[5,3,2]; % number of inputs in each stage
    A=[1.2,1.3,1.1];
    B=[1.5,1.1,0.8];
    gamma=[0.08, -0.20, 0.46];
    theta=[0.58,0.67,0.73];
    
    alpha1=[0.25,0.16,0.04,0.35,0.20];
    alpha2=[0.35,0.35,0.30];
    alpha3=[0.8,0.2];
    
    ratio0=[1.0,1.0,1.0]; %% combination coefficients of y0
    ratio1=[0.5,0.5;... %% combination coefficients of realized stage 0 weather
        0.4,0.6;... %% expect period 1 weather, combination coefficients of realized stage 0 weather
        0.3,0.7;... %% expect period 2 weather, combination coefficients of realized stage 0 weather
        0.2,0.8];   %% expect period 3 weather, combination coefficients of realized stage 0 weather
    ratio2=[0.5,0.5;... %% combination coefficients of realized stage 1 weather
        0.1,0.1;0.3,0.7;... %% expect period 2 weather, combination coefficients of realized stage 0 weather, stage 1 weather
        0.1,0.1;0.2,0.8];   %% expect period 3 weather, combination coefficients of realized stage 0 weather, stage 1 weather
    
    ratio3=[0.5,0.5;... %% combination coefficients of realized stage 2 weather
        0.0,0.0;0.1,0.1;0.2,0.8]; %% expect period 3 weather, combination coefficients of realized stage 0 weather, stage 1 weather, stage 2 weather

    
end

