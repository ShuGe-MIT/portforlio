function [supx, probx] = ini_normal(meanx,stdx,n)
    % this function initialize a discretized normal distribution given the
    % mean and std of a normal distribution
    minx=meanx-3*stdx;
    maxx=meanx+3*stdx;
    minx=max(minx,0);
    
    gridx=linspace(minx,maxx,(n+1));
    cdfx=normcdf(gridx,meanx,stdx);
    
    probx=diff(cdfx);
    probx=probx./sum(probx);
    supx(1:n)=0.0;
    for i=1:n
        supx(i)=(gridx(i)+gridx(i+1))/2;
    end
end

