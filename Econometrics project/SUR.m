function [BETA, COV_inv]=SUR(no_equations, obs, log_initial,  log_price, log_raintemp, lnxy_dist_cal_sai)
    % this function does seemingly unrelated regression according to
    % section 2 of SUR_robust_cluster_vce_2020_05_25.pdf.
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
%     SUR_load_data=toc

    %% OLS
    tic;
    eps_hat_by_eq(1:obs, 1:no_equations)=0.0;
    eps_hat_by_eq(1:obs, 1)=Y_ori(:,1)-X_ori1*((X_ori1'*X_ori1)\X_ori1'*Y_ori(:,1));
    eps_hat_by_eq(1:obs, 2:(N(3)+1))=Y_ori(:,2:(N(3)+1))-X_ori2*((X_ori2'*X_ori2)\X_ori2'*Y_ori(:,2:(N(3)+1)));
    eps_hat_by_eq(1:obs, (N(3)+2):(N(3)+N(2)+1))=Y_ori(:,(N(3)+2):(N(3)+N(2)+1))-X_ori3*((X_ori3'*X_ori3)\X_ori3'*Y_ori(:,(N(3)+2):(N(3)+N(2)+1)));
    eps_hat_by_eq(1:obs, (N(3)+N(2)+2):(N(3)+N(2)+N(1)+1))=Y_ori(:,(N(3)+N(2)+2):(N(3)+N(2)+N(1)+1))-X_ori4*((X_ori4'*X_ori4)\X_ori4'*Y_ori(:,(N(3)+N(2)+2):(N(3)+N(2)+N(1)+1)));

    cov=(eps_hat_by_eq')*eps_hat_by_eq./obs;
    cov_inv=inv(cov);
%     SUR_OLS=toc    
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
    BETA=COV_inv\(temp1*Y);
%     SUR_GLS=toc
end

