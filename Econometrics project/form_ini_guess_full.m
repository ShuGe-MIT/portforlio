for www=1:5
    clear;
%     load('result_with_regression_little_sai.mat','A','B','ratio0','ratio1','ratio2','ratio3');
%     load('guess_initial_alpha.mat');
    MC=512;
    guess(1:MC,1:51)=0.0;
    %% First, let's form the guess for A
    for i=1:MC
        guess(i,17:19)=random_num(1,3,0.1,5.0);
    end

    %% Second, let's form the guess for B
    for i=1:MC
        guess(i,20:22)=random_num(1,3,0.1,5.0);
    end

    %% Third, let's form the guess for ratio and y0
    for i=1:MC
        guess(i,23:51)=random_num(1,29,0,5);
%         guess(i,23:51)=[ratio1(1,:), ratio1(2,:), ratio1(3,:), ratio1(4,:),...
%         ratio2(1,:), ratio2(2,:), ratio2(3,:), ratio2(4,:), ratio2(5,:),...
%         ratio3(1,:), ratio3(2,:), ratio3(3,:), ratio3(4,:), initial_guess]; 
    end

    %% Fourth, let's form the guess for alpha
    for i=1:MC
        guess(i,15:16)=rand(1,2); %% alpha3
        guess(i,12:14)=rand(1,3); %% alpha2
        guess(i,7:11)=rand(1,5); %% alpha1
    end

    %% Lastly, let's form the guess for gamma, theta
    guess_a(1:512,1:6)=0.0;
    interval_gamma=[-10,-5;-5,-1;-1,0;0,1.0];
    interval_theta=[0,0.5;0.5,1];
    for i=1:511
        signal=de2bi(i,9);
        signal_gamma1=bi2de(signal(1:2))+1;
        signal_gamma2=bi2de(signal(3:4))+1;
        signal_gamma3=bi2de(signal(5:6))+1;
        signal_theta1=bi2de(signal(7))+1;
        signal_theta2=bi2de(signal(8))+1;
        signal_theta3=bi2de(signal(9))+1;
        guess_a(i,1)=random_num(1,1,interval_gamma(signal_gamma1,1),interval_gamma(signal_gamma1,2)); % gamma1
        guess_a(i,2)=random_num(1,1,interval_gamma(signal_gamma2,1),interval_gamma(signal_gamma2,2)); % gamma2
        guess_a(i,3)=random_num(1,1,interval_gamma(signal_gamma3,1),interval_gamma(signal_gamma3,2)); % gamma3
        guess_a(i,4)=random_num(1,1,interval_theta(signal_theta1,1),interval_theta(signal_theta1,2)); % theta1
        guess_a(i,5)=random_num(1,1,interval_theta(signal_theta2,1),interval_theta(signal_theta2,2)); % theta2
        guess_a(i,6)=random_num(1,1,interval_theta(signal_theta3,1),interval_theta(signal_theta3,2)); % theta3
    end
    i=512;
    signal=zeros(1,9);
    signal_gamma1=bi2de(signal(1:2))+1;
    signal_gamma2=bi2de(signal(3:4))+1;
    signal_gamma3=bi2de(signal(5:6))+1;
    signal_theta1=bi2de(signal(7))+1;
    signal_theta2=bi2de(signal(8))+1;
    signal_theta3=bi2de(signal(9))+1;
    guess_a(i,1)=random_num(1,1,interval_gamma(signal_gamma1,1),interval_gamma(signal_gamma1,2)); % gamma1
    guess_a(i,2)=random_num(1,1,interval_gamma(signal_gamma2,1),interval_gamma(signal_gamma2,2)); % gamma2
    guess_a(i,3)=random_num(1,1,interval_gamma(signal_gamma3,1),interval_gamma(signal_gamma3,2)); % gamma3
    guess_a(i,4)=random_num(1,1,interval_theta(signal_theta1,1),interval_theta(signal_theta1,2)); % theta1
    guess_a(i,5)=random_num(1,1,interval_theta(signal_theta2,1),interval_theta(signal_theta2,2)); % theta2
    guess_a(i,6)=random_num(1,1,interval_theta(signal_theta3,1),interval_theta(signal_theta3,2)); % theta3

    %% To test that my algorithm is correct
%     sum(guess_a(:,1)<-5);
%     sum(guess_a(:,1)>-5);
%     sum(guess_a(:,2)<-5);
%     sum(guess_a(:,2)>-5);
%     sum(guess_a(:,3)<-5);
%     sum(guess_a(:,3)>-5);
%     sum(guess_a(:,4)<0.5);
%     sum(guess_a(:,4)>0.5);
%     sum(guess_a(:,5)<0.5);
%     sum(guess_a(:,5)>0.5);
%     sum(guess_a(:,6)<0.5);
%     sum(guess_a(:,6)>0.5);

    %% Finally, form the whole matrix
    guess(:,1:6)=guess_a(:,1:6);
    % save('guess.mat');
    kk=zeros(2000,13);
    
    load('data.mat');
    valid(1:MC)=0;
    for i=1:MC
        [xy_dist_cal]=solve_x(guess(i,:), data(:,1:30));
        if xy_dist_cal==real(xy_dist_cal)
            if isinf(xy_dist_cal)==kk
                valid(i)=1;
            end
        end
    end
    guess(logical(1-valid),:) = [];
    save('guess_valid_new.mat','guess');
    a=guess;
    load('guess_valid_new.mat', 'guess');
    guess=[a;guess];
    save('guess_valid_new.mat','guess');

end