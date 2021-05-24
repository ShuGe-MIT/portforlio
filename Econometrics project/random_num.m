function [result] = random_num(x,y,lb,ub)
    %%% This function generates random number from uniform dist matrix size x*y with lb and
    %%% ub
    result(1:x,1:y)=rand(x,y);
    for i=1:x
        for j=1:y
            result(i,j)=(ub-lb)*result(i,j)+lb;
        end
    end
end

