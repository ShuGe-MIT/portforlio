function [BETA, COV_inv]=SUR(no_equations, obs, log_initial,  log_price, log_raintemp, lnxy_dist_cal_sai)
tic;    
warning('off');
    [~, obs, n_cluster, cum_cluster, count_cluster, ~, N, ~, ~,...
    ~, ~,...
    ~, ~, ~, ~, ~, ~, ~]=parameter();
    
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
time_for_SUR_matrix=toc

end

