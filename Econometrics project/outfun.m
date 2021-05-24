

function stop = outfun(x, optimValues, state)
% this function sets the stopping condition for fmincon
stop = (optimValues.fval<300);
end

