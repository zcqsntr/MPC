function c = cost(xk, e, w, y)

    target_N = [250, 550];
    %COST Summary of this function goes here
    %   Detailed explanation goes here
    c = sum(sum(((target_N-xk(2:end,1:2))/1000).^2));
    
end

