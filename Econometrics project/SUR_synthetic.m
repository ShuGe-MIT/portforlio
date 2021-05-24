function [BETA, COV]=SUR_synthetic(no_equations, obs, log_initial,  log_price, log_raintemp, lnxy_dist_cal_sai)
    % this function does seemingly unrelated regression according to
    % section 3.2 of SUR_robust_cluster_vce_2020_05_25.pdf.
    %tic;    
    warning('off');
    [~, obs, n_cluster, cum_cluster, count_cluster, ~, N, ~, ~,...
    ~, ~,...
    ~, ~, ~, ~, ~, ~, ~]=parameter();

%     H=max(size(count_cluster)); % number of strata
%     Q=max(size(n_cluster)); % number of clusters in total
%     Q_h=count_cluster; % number of clusters in each stratum
%     P_qh=n_cluster; % number of observations in each cluster
%     cum_strata=cumsum(count_cluster);
    
    
    variables=[length(log_initial(1,:))+length(log_price(1,:))-1+length(log_raintemp(1,:))+1];
    no_var=[variables,ones(1,N(3))*(variables-2),ones(1,N(2))*(variables-3),ones(1,N(1))*(variables-5)];
    no_var_cum=cumsum(no_var);
    no_var_cum=[0, no_var_cum];
    K=sum(no_var);
    
    %% data processing
    Y_ori(1:obs,1:no_equations)=[lnxy_dist_cal_sai(:,(N(1)+N(2)+N(3))+3), lnxy_dist_cal_sai(:,(N(1)+N(2)+1):(N(1)+N(2)+N(3))), ...
                                lnxy_dist_cal_sai(:,(N(1)+1):(N(1)+N(2))), lnxy_dist_cal_sai(:,1:N(1))];
    a=ones(obs, 1);
%     save('X_log_initial','log_initial');
    lnxy_dist_cal_sai_=lnxy_dist_cal_sai(:,1:(N(1)+N(2)+N(3)));
%     save('X_lnxy_dist_cal_sai_','lnxy_dist_cal_sai_');
    X_ori1=[log_initial, lnxy_dist_cal_sai(:,1:(N(1)+N(2)+N(3))), log_raintemp, a];
    X_ori2=[log_initial, lnxy_dist_cal_sai(:,1:(N(1)+N(2))), log_price(:,(N(1)+N(2)+1):(N(1)+N(2)+N(3)))-log_price(:,(N(1)+N(2)+N(3))+1), log_raintemp(:,1:6), a];
    X_ori3=[log_initial, lnxy_dist_cal_sai(:,1:N(1)), log_price(:,(N(1)+N(2)+1):(N(1)+N(2)+N(3))), log_price(:,(N(1)+1):(N(1)+N(2))), log_price(:,(N(1)+N(2)+N(3))+1), log_raintemp(:,1:4), a];
    X_ori4=[log_initial, log_price(:,(N(1)+N(2)+1):(N(1)+N(2)+N(3))), log_price(:,(N(1)+1):(N(1)+N(2))), log_price(:,1:N(1)), log_price(:,(N(1)+N(2)+N(3))+1), log_raintemp(:,1:2), a];
    
    Y=Y_ori(:);
    X=blkdiag(X_ori1, X_ori2, X_ori2, X_ori3, X_ori3, X_ori3, X_ori4, X_ori4, X_ori4, X_ori4, X_ori4);
%     save('X','X');
%     save('Y','Y');
    % SUR_syn_load=toc

    %% OLS
    % tic;
    eps_hat_by_eq(1:obs, 1:no_equations)=0.0;   
    Y_oriii=Y_ori(:,1:3);
%     save('Y_ori1','Y_oriii');
%     save('X_ori1','X_ori1');
%     save('X_ori2','X_ori2');
    eps_hat_by_eq(1:obs, 1)=Y_ori(:,1)-X_ori1*((X_ori1'*X_ori1)\X_ori1'*Y_ori(:,1));
    eps_hat_by_eq(1:obs, 2:(N(3)+1))=Y_ori(:,2:(N(3)+1))-X_ori2*((X_ori2'*X_ori2)\X_ori2'*Y_ori(:,2:(N(3)+1)));
    eps_hat_by_eq(1:obs, (N(3)+2):(N(3)+N(2)+1))=Y_ori(:,(N(3)+2):(N(3)+N(2)+1))-X_ori3*((X_ori3'*X_ori3)\X_ori3'*Y_ori(:,(N(3)+2):(N(3)+N(2)+1)));
    eps_hat_by_eq(1:obs, (N(3)+N(2)+2):(N(3)+N(2)+N(1)+1))=Y_ori(:,(N(3)+N(2)+2):(N(3)+N(2)+N(1)+1))-X_ori4*((X_ori4'*X_ori4)\X_ori4'*Y_ori(:,(N(3)+N(2)+2):(N(3)+N(2)+N(1)+1)));
    cov=(eps_hat_by_eq')*eps_hat_by_eq./obs;
%     save('eps_hat_by_eq','eps_hat_by_eq');
    cov_inv=inv(cov);
    % SUR_syn_OLS=toc    
    %% GLS
%     tic;    
    temp1(1:sum(no_var),1:(no_equations*obs))=0.0;
    
    for i=1:no_equations
        for j=1:no_equations
            temp1((no_var_cum(i)+1):(no_var_cum(i+1)),((j-1)*obs+1):(j*obs))=...
                X(((i-1)*obs+1):(i*obs),(no_var_cum(i)+1):(no_var_cum(i+1)))'.*cov_inv(i,j);
        end
    end
    COV_inv=temp1*X;
    BETA=COV_inv\(temp1*Y);
    
%     SUR_syn_GLS=toc
    %% GLS
%     R=no_equations;
%     K=no_var_cum(R+1);
%     X(1:(R*obs),1:K)=0.0;
%     for i=1:obs
%         X((i-1)*R+1:i*R,:)=blkdiag(X_ori1(i,:), X_ori2(i,:), X_ori2(i,:),...
%             X_ori3(i,:), X_ori3(i,:), X_ori3(i,:), X_ori4(i,:), X_ori4(i,:), X_ori4(i,:), X_ori4(i,:), X_ori4(i,:));
%     end
%             
%     Y=Y_ori';
%     Y=Y(:);
%     temp1(1:K,1:R*obs)=X'*(kron(eye(obs),cov_inv));
%     COV_inv=temp1*X;
%     BETA=COV_inv\(temp1*Y);
    %% Robust covariance without clustering
    % tic;   
    e_check=Y-X*BETA;
    e_check=reshape(e_check,[obs,no_equations]);
    Omega_check(1:no_equations,1:no_equations)=e_check'*e_check./obs;
    Omega_check_inv=inv(Omega_check);
    
    
    D_hat(1:K,1:K)=0.0;
    GG_hat(1:K,1:K)=0.0;
    temp2=0.0;
    temp3=0.0;
%     tic;
    for i=1:obs
        %x_n=blkdiag(X_ori1(i,:), X_ori2(i,:), X_ori2(i,:), X_ori3(i,:), X_ori3(i,:), X_ori3(i,:),...
        %    X_ori4(i,:), X_ori4(i,:), X_ori4(i,:), X_ori4(i,:), X_ori4(i,:));
        help=i:obs:(i+obs*(no_equations-1));
        x_n=X(help,:);
        temp2=x_n'*Omega_check_inv;
%          save('Omega_check_inv','Omega_check_inv')
%         temp2(1:K,1:no_equations)=0.0;
%         for j=1:no_equations
%             for k=1:no_equations
%                 temp2((no_var_cum(j)+1):(no_var_cum(j+1)),k)=...
%                     x_n(j,(no_var_cum(j)+1):(no_var_cum(j+1)))'.*Omega_check_inv(j,k);
%             end
%         end
        temp3=temp2*e_check(i,:)';
%        save('e_check','e_check');
        D_hat=D_hat+temp2*x_n;
%        save('temp2','temp2');
%        save('x_n','x_n');
        GG_hat=GG_hat+temp3*(temp3)';
%        save('temp3','temp3');
        
    end
    % toc;
    D_hat=inv(D_hat);
%    save('D_hat','D_hat');
    GG_hat=GG_hat.*(obs/(obs-K));
%    save('GG_hat','GG_hat');
    COV=D_hat*GG_hat*D_hat;
%    save('COV','COV');
%     eig(COV)
%     robust_covariance=toc
    
    %% Allowing for clustering
%     R=no_equations;
%     Z(1:(R*obs),1:obs)=0.0;
%     for i=1:obs
%         Z(((i-1)*R+1):i*R,i)=1;
%     end
%     I_RRN=kron(ones(1,obs),eye(R));
%     e_check=Y-X*BETA;
%     eps_check=kron(ones(1,obs),e_check).*Z;
%     OMEGA_check(1:no_equations,1:no_equations)=(I_RRN*(eps_check*eps_check')*I_RNR)./obs;
%     OMEGA_check_inv=inv(OMEGA_check);
%     
%     GG_hat(1:K,1:K)=0.0;
%     D_hat=inv(X'*(kron(eye(obs),OMEGA_check_inv))*X);
%     
%     GG_hat=X'*kron(eye(obs),OMEGA_check_inv)*eps_check;
%     GG_hat=GG_hat*GG_hat';
%     GG_hat=GG_hat.*(obs/(obs-K));
%     COV=D_hat*GG_hat*D_hat;

end

