function [DELTA, COV]=SUR(no_equations, obs, log_initial,  log_price, log_raintemp, lnxy_dist_cal_sai)
    % this function does seemingly unrelated regression according to
    % section 3.1 of SUR_robust_cluster_vce_2020_05_25.pdf.
    tic;
    warning('off');
    [~, obs, n_cluster, cum_cluster, count_cluster, ~, N, ~, ~,...
    ~, ~,...
    ~, ~, ~, ~, ~, ~, ~]=parameter();
    
    H=max(size(count_cluster)); % number of strata
    Q=max(size(n_cluster)); % number of clusters in total
    Q_h=count_cluster; % number of clusters in each stratum
    P_qh=n_cluster; % number of observations in each cluster
    cum_strata=cumsum(count_cluster);
    
    
    variables=[length(log_initial(1,:))+length(log_price(1,:))-1+length(log_raintemp(1,:))+1];
    no_var=[variables,ones(1,N(3))*(variables-2),ones(1,N(2))*(variables-3),ones(1,N(1))*(variables-5)];
    no_var_cum=cumsum(no_var);
    no_var_cum=[0, no_var_cum];
    K=sum(no_var);
    
    %% data processing
    Y_ori(1:obs,1:no_equations)=[lnxy_dist_cal_sai(:,(N(1)+N(2)+N(3))+3), lnxy_dist_cal_sai(:,(N(1)+N(2)+1):(N(1)+N(2)+N(3))), ...
                                lnxy_dist_cal_sai(:,(N(1)+1):(N(1)+N(2))), lnxy_dist_cal_sai(:,1:N(1))];
    a=ones(obs, 1);
    
    X_ori1=[log_initial, lnxy_dist_cal_sai(:,1:(N(1)+N(2)+N(3))), log_raintemp, a];
    X_ori2=[log_initial, lnxy_dist_cal_sai(:,1:(N(1)+N(2))), log_price(:,(N(1)+N(2)+1):(N(1)+N(2)+N(3)))-log_price(:,(N(1)+N(2)+N(3))+1), log_raintemp(:,1:6), a];
    X_ori3=[log_initial, lnxy_dist_cal_sai(:,1:N(1)), log_price(:,(N(1)+N(2)+1):(N(1)+N(2)+N(3))), log_price(:,(N(1)+1):(N(1)+N(2))), log_price(:,(N(1)+N(2)+N(3))+1), log_raintemp(:,1:4), a];
    X_ori4=[log_initial, log_price(:,(N(1)+N(2)+1):(N(1)+N(2)+N(3))), log_price(:,(N(1)+1):(N(1)+N(2))), log_price(:,1:N(1)), log_price(:,(N(1)+N(2)+N(3))+1), log_raintemp(:,1:2), a];
    
    Y=Y_ori(:);
    X=blkdiag(X_ori1, X_ori2, X_ori2, X_ori3, X_ori3, X_ori3, X_ori4, X_ori4, X_ori4, X_ori4, X_ori4);
%     SUR_real_loading=toc
    %% OLS
    tic;    
    eps_hat_by_eq(1:obs, 1:no_equations)=0.0;
    eps_hat_by_eq(1:obs, 1)=Y_ori(:,1)-X_ori1*((X_ori1'*X_ori1)\X_ori1'*Y_ori(:,1));
    eps_hat_by_eq(1:obs, 2:(N(3)+1))=Y_ori(:,2:(N(3)+1))-X_ori2*((X_ori2'*X_ori2)\X_ori2'*Y_ori(:,2:(N(3)+1)));
    eps_hat_by_eq(1:obs, (N(3)+2):(N(3)+N(2)+1))=Y_ori(:,(N(3)+2):(N(3)+N(2)+1))-X_ori3*((X_ori3'*X_ori3)\X_ori3'*Y_ori(:,(N(3)+2):(N(3)+N(2)+1)));
    eps_hat_by_eq(1:obs, (N(3)+N(2)+2):(N(3)+N(2)+N(1)+1))=Y_ori(:,(N(3)+N(2)+2):(N(3)+N(2)+N(1)+1))-X_ori4*((X_ori4'*X_ori4)\X_ori4'*Y_ori(:,(N(3)+N(2)+2):(N(3)+N(2)+N(1)+1)));

    cov=(eps_hat_by_eq')*eps_hat_by_eq./obs;
    cov_inv=inv(cov);
%     SUR_real_OLS=toc    
    %% GLS
    tic;    
    temp1(1:sum(no_var),1:(no_equations*obs))=0.0;
    
    for i=1:no_equations
        for j=1:no_equations
            temp1((no_var_cum(i)+1):(no_var_cum(i+1)),((j-1)*obs+1):(j*obs))=...
                X(((i-1)*obs+1):(i*obs),(no_var_cum(i)+1):(no_var_cum(i+1)))'.*cov_inv(i,j);
        end
    end
    COV_inv=temp1*X;
    DELTA=COV_inv\(temp1*Y);
%     SUR_real_GLS=toc    
    %% Allowing for clustering
    tic;    
    e_check=Y-X*DELTA;
    e_check=reshape(e_check,[obs,no_equations]);
    Omega_check(1:no_equations,1:no_equations)=e_check'*e_check./obs;
    Omega_check_inv=inv(Omega_check);
%     time1=toc
    tic;
    u_pqh(1:obs,1:K)=0.0;
    for i=1:obs
        x_pqh=blkdiag(X_ori1(i,:), X_ori2(i,:), X_ori2(i,:), X_ori3(i,:), X_ori3(i,:), X_ori3(i,:),...
            X_ori4(i,:), X_ori4(i,:), X_ori4(i,:), X_ori4(i,:), X_ori4(i,:));
        u_pqh(i,:)=(x_pqh')*(Omega_check_inv)*(e_check(i,:)');
    end
%     SUR_real_u_pqh=toc
    tic;
    u_qh(1:Q, 1:K)=0.0;
    u_qh(1,:)=sum(u_pqh(1:cum_cluster(1),:));
    for i=2:Q
        u_qh(i,:)=sum(u_pqh((cum_cluster(i-1)+1):cum_cluster(i),:));
    end
%     SUR_real_u_qh=toc
    tic;
    u_h(1:H,1:K)=0.0;
    u_h(1,:)=mean(u_qh(1:cum_strata(1),:));
    for i=2:H
        u_h(i,:)=mean(u_qh((cum_strata(i-1)+1):(cum_strata(i)),:));
    end
%     SUR_real_u_h=toc
    tic;
    G_hat(1:K,1:K)=0.0;
    temp=0.0;
    no=1;
%     udiff(1:Q,1:K)=0.0;
%     udiff2(1:Q,1:K,1:K)=0.0;
%     bracket(1:H,1:K,1:K)=0.0;
    for i=1:H
        for j=1:Q_h(i)
            temp=temp+((u_qh(no,:)-u_h(i,:))')*(u_qh(no,:)-u_h(i,:));
%             udiff(no,:)=((u_qh(no,:)-u_h(i,:))');
%             udiff2(no,:,:)=((u_qh(no,:)-u_h(i,:))')*(u_qh(no,:)-u_h(i,:));
            no=no+1;
        end
%         bracket(i,:,:)=temp;
        G_hat=G_hat+temp.*(Q_h(i)/(Q_h(i)-1));
        temp=0.0;
    end
%     save('udiff.mat','H','Q','K','Q_h','P_qh','udiff','udiff2','bracket');
    G_hat=G_hat.*((obs-1)/(obs-K));
%     SUR_real_G_hat=toc
    tic;
    D_hat(1:K,1:K)=0.0;
    temp2(1:sum(no_var),1:(no_equations*obs))=0.0;
    
    for i=1:no_equations
        for j=1:no_equations
            temp2((no_var_cum(i)+1):(no_var_cum(i+1)),((j-1)*obs+1):(j*obs))=...
                X(((i-1)*obs+1):(i*obs),(no_var_cum(i)+1):(no_var_cum(i+1)))'.*Omega_check_inv(i,j);
        end
    end
    D_hat=inv(temp2*X);
%     SUR_real_D_hat=toc
    tic;
    %D_hat=inv(X'*(kron(Omega_check_inv,eye(obs)))*X);
    COV=D_hat*G_hat*D_hat;
%     SUR_real_cov=toc
end

