function [BETA, COV]=SUR(no_equations, obs, log_initial,  log_price, log_raintemp, lnxy_dist_cal_sai)
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

    %% OLS
    BETA1=(X_ori1'*X_ori1)\X_ori1'*Y_ori(:,1);
    BETA2=(X_ori2'*X_ori2)\X_ori2'*Y_ori(:,2:(N(3)+1));
    BETA3=(X_ori3'*X_ori3)\X_ori3'*Y_ori(:,(N(3)+2):(N(3)+N(2)+1));
    BETA4=(X_ori4'*X_ori4)\X_ori4'*Y_ori(:,(N(3)+N(2)+2):(N(3)+N(2)+N(1)+1));
    DELTA_OLS=[BETA1;BETA2(:);BETA3(:);BETA4(:)];
    e=Y-X*DELTA_OLS;
    e=reshape(e,[obs,no_equations]);
    e=e';
    e=e(:);
    R=no_equations;
    Z(1:(R*obs),1:obs)=0.0;
    for i=1:obs
        Z(((i-1)*R+1):i*R,i)=1;
    end
    
    eps=kron(ones(1,obs),e).*Z; %% kron(ones(1,obs),e)=kron(e,ones(1,obs))
    I_RRN=kron(ones(1,obs),eye(R));
    I_RNR=kron(ones(obs,1),eye(R));
    cov=(I_RRN*(eps*eps')*I_RNR)./obs;
    cov_inv=inv(cov);
    
    %% GLS
    K=no_var_cum(R+1);
    x(1:(R*obs),1:K)=0.0;
    for i=1:obs
        x((i-1)*R+1:i*R,:)=blkdiag(X_ori1(i,:), X_ori2(i,:), X_ori2(i,:),...
            X_ori3(i,:), X_ori3(i,:), X_ori3(i,:), X_ori4(i,:), X_ori4(i,:), X_ori4(i,:), X_ori4(i,:), X_ori4(i,:));
    end
            
    y=Y_ori';
    y=y(:);
    temp1(1:K,1:R*obs)=x'*(kron(eye(obs),cov_inv));
    COV_inv=temp1*x;
    BETA=COV_inv\(temp1*y);
    
    %% Allowing for clustering
    e_check=y-x*BETA;
    eps_check=kron(ones(1,obs),e_check).*Z;
    OMEGA_check(1:no_equations,1:no_equations)=(I_RRN*(eps_check*eps_check')*I_RNR)./obs;
    OMEGA_check_inv=inv(OMEGA_check);
    
%     u_bqh(1:obs,1:K)=0.0;
%     e_check_temp=reshape(e_check,[R, obs]);
%     for i=1:obs
%         x_bqh=blkdiag(X_ori1(i,:), X_ori2(i,:), X_ori2(i,:), X_ori3(i,:), X_ori3(i,:), X_ori3(i,:),...
%             X_ori4(i,:), X_ori4(i,:), X_ori4(i,:), X_ori4(i,:), X_ori4(i,:));
%         u_bqh(i,:)=(x_bqh')*(OMEGA_check_inv)*(e_check_temp(:,i));
%     end
    
    
    u_qh(1:Q, 1:K)=0.0;
    u_qh(1,:)=x(1:cum_cluster(1)*R,:)'*kron(eye(cum_cluster(1)),OMEGA_check_inv)*e_check(1:cum_cluster(1)*R);
    for i=2:Q
        u_qh(i,:)=x((cum_cluster(i-1)*R+1):cum_cluster(i)*R,:)'*kron(eye(n_cluster(i)),OMEGA_check_inv)*e_check((cum_cluster(i-1)*R+1):cum_cluster(i)*R);
    end
    
    u_h_bar(1:H,1:K)=0.0;
    u_h_bar(1,:)=mean(u_qh(1:cum_strata(1),:));
    for i=2:H
        u_h_bar(i,:)=mean(u_qh((cum_strata(i-1)+1):(cum_strata(i)),:));
    end
    
    G_hat(1:K,1:K)=0.0;
    u_h_tilda=u_qh(1:cum_strata(1),:)';
    u_h_tilda=u_h_tilda(:);
    u_h_bracket=u_h_tilda-kron(ones(Q_h(1),1),u_h_bar(1,:)');
    Z_QhKQh=0;
    Z_QhKQh(1:Q_h(1)*K,Q_h(1))=0.0;
    for i=1:Q_h(1)
        Z_QhKQh(((i-1)*K+1):i*K,i)=1;
    end
    U_h=kron(ones(1,Q_h(1)),u_h_bracket).*Z_QhKQh;
    I_KQhK=kron(ones(1,Q_h(1)),eye(K));
    I_QhKK=kron(ones(Q_h(1),1),eye(K));
    bracket=I_KQhK*(U_h*U_h')*I_QhKK;
    G_hat=G_hat+Q_h(1)/(Q_h(1)-1)*bracket;
    for i=2:H
        u_h_tilda=u_qh(cum_strata(i-1)+1:cum_strata(i),:)';
        u_h_tilda=u_h_tilda(:);
        u_h_bracket=u_h_tilda-kron(ones(Q_h(i),1),u_h_bar(i,:)');
        Z_QhKQh=0;
        Z_QhKQh(1:Q_h(i)*K,Q_h(i))=0.0;
        for j=1:Q_h(i)
            Z_QhKQh(((j-1)*K+1):j*K,j)=1;
        end
        U_h=kron(ones(1,Q_h(i)),u_h_bracket).*Z_QhKQh;
        I_KQhK=kron(ones(1,Q_h(i)),eye(K));
        I_QhKK=kron(ones(Q_h(i),1),eye(K));
        bracket=I_KQhK*(U_h*U_h')*I_QhKK;
        G_hat=G_hat+Q_h(i)/(Q_h(i)-1)*bracket;
    end
    G_hat=G_hat.*((obs-1)/(obs-K));
   
    D_hat=inv(x'*(kron(eye(obs),OMEGA_check_inv))*x);
    COV=D_hat*G_hat*D_hat;
end

