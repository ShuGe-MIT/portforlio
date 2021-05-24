% This file polishes the initial guess obtained from form_ini_guess_full.m.
% For each set of structural parameter values, it tries to update it such
% that the Wald metric is closer to 0.

clear;
load('guess_valid_new.mat');
load('data.mat');
load('result_with_regression_little_sai.mat','sai_dist','lnxy_dist_sai');
no_equations=11;
obs=2000;
log_initial=data(:,50:52);
log_price_dist=data(:,31:41);
log_rain_temp_dist=data(:,42:49);

[BETA_data, W_data]=SUR(no_equations, obs, log_initial,  log_price_dist, log_rain_temp_dist, lnxy_dist_sai);

test_start=1;
test_end=200;
guess_polished(test_start:test_end,1:51)=0.0;
for i=test_start:test_end
    i    
    tic;
    para_guess=guess(i,:);
    temp=squeeze(sai_dist(i,:,:,:));
    
    guess_candidate(1:64,1:51)=0.0;
    error(1:64)=0.0;
    for j=1:64
        guess_candidate(j,1:51)=para_guess;
        signal=de2bi(j-1,6);
        if signal(1)==1
            guess_candidate(j,1)=-guess_candidate(j,1);
        end
        if signal(2)==1
            guess_candidate(j,2)=-guess_candidate(j,2);
        end
        if signal(3)==1
            guess_candidate(j,3)=-guess_candidate(j,3);
        end
        if signal(4)==1
            guess_candidate(j,4)=0.5;
        end
        if signal(5)==1
            guess_candidate(j,5)=0.5;
        end
        if signal(6)==1
            guess_candidate(j,6)=0.5;
        end
        error(j)=compute_error(guess_candidate(j,:), temp, data, BETA_data, W_data);
        if (error(j)~=real(error(j)))
            error(j)=1e10;    
        end
    end
    
    [~,help]=min(error);
    guess_polished(i,:)=guess_candidate(help,:);
    toc;
end

save('guess_polished.mat','guess_polished_11_22');