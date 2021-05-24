function [w1_sup, w2_sup, w3_sup, p_sup, w1_prob, w2_prob, w3_prob, p_prob] = price(n)
    % this function constructs the distribution of the prices

    [w1_sup(1,1,:), w1_prob(1,1,:)]=ini_normal(10.71, 4.44, n); % price of chemical in first quarter
%     calculate_std(squeeze(w1_sup(1,1,:)),squeeze(w1_prob(1,1,:)),n);
    [w1_sup(1,2,:), w1_prob(1,2,:)]=ini_normal(8.09, 3.64, n); % price of seeds in first quarter
%     calculate_std(squeeze(w1_sup(1,2,:)),squeeze(w1_prob(1,2,:)),n);
    [w1_sup(1,3,:), w1_prob(1,3,:)]=ini_normal(1.70, 2.66, n); % price of seedlings in first quarter
%     calculate_std(squeeze(w1_sup(1,3,:)),squeeze(w1_prob(1,3,:)),n);
    [w1_sup(1,4,:), w1_prob(1,4,:)]=ini_normal(101.37, 123.39, n); % price of labor in first quarter
%     calculate_std(squeeze(w1_sup(1,4,:)),squeeze(w1_prob(1,4,:)),n);
    [w1_sup(1,5,:), w1_prob(1,5,:)]=ini_normal(120.0, 23.39, n); % price of equipment in first quarter
%     calculate_std(squeeze(w1_sup(1,5,:)),squeeze(w1_prob(1,5,:)),n);

    w1_sup(2,:,:)=w1_sup(1,:,:);
    w1_sup(3,:,:)=w1_sup(1,:,:);
    w1_sup(4,:,:)=w1_sup(1,:,:);
    w1_prob(2,:,:)=w1_prob(1,:,:);
    w1_prob(3,:,:)=w1_prob(1,:,:);
    w1_prob(4,:,:)=w1_prob(1,:,:);
    
    [w2_sup(1,1,:), w2_prob(1,1,:)]=ini_normal(10.76, 4.42, n); % price of chemical in first quarter
%     calculate_std(squeeze(w2_sup(1,1,:)),squeeze(w2_prob(1,1,:)),n);
    [w2_sup(1,2,:), w2_prob(1,2,:)]=ini_normal(69.06, 67.26, n); % price of labor in first quarter
%     calculate_std(squeeze(w2_sup(1,2,:)),squeeze(w2_prob(1,2,:)),n);
    [w2_sup(1,3,:), w2_prob(1,3,:)]=ini_normal(80.0, 27.26, n); % price of equipment in first quarter
%     calculate_std(squeeze(w2_sup(1,3,:)),squeeze(w2_prob(1,3,:)),n);

    w2_sup(2,:,:)=w2_sup(1,:,:);
    w2_sup(3,:,:)=w2_sup(1,:,:);
    w2_sup(4,:,:)=w2_sup(1,:,:);
    w2_prob(2,:,:)=w2_prob(1,:,:);
    w2_prob(3,:,:)=w2_prob(1,:,:);
    w2_prob(4,:,:)=w2_prob(1,:,:);
    
    [w3_sup(1,1,:), w3_prob(1,1,:)]=ini_normal(127.61, 95.79, n); % price of labor in first quarter
%     calculate_std(squeeze(w3_sup(1,1,:)),squeeze(w3_prob(1,1.:)),n);
    [w3_sup(1,2,:), w3_prob(1,2,:)]=ini_normal(140.0, 55.79, n); % price of equipment in first quarter
%     calculate_std(squeeze(w3_sup(1,2,:)),squeeze(w3_prob(1,2,:)),n);

    w3_sup(2,:,:)=w3_sup(1,:,:);
    w3_sup(3,:,:)=w3_sup(1,:,:);
    w3_sup(4,:,:)=w3_sup(1,:,:);
    w3_prob(2,:,:)=w3_prob(1,:,:);
    w3_prob(3,:,:)=w3_prob(1,:,:);
    w3_prob(4,:,:)=w3_prob(1,:,:);
    
    [p_sup(1,:), p_prob(1,:)]=ini_normal(7.35, 2.73, n); % price of rice
%     calculate_std(squeeze(p_sup(1,:)),squeeze(p_prob(1,:)),n);
    p_sup(2,:)=p_sup(1,:);
    p_sup(3,:)=p_sup(1,:);
    p_sup(4,:)=p_sup(1,:);
    p_prob(2,:)=p_prob(1,:);
    p_prob(3,:)=p_prob(1,:);
    p_prob(4,:)=p_prob(1,:);
    

end

