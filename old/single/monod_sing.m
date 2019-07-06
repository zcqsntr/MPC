function growth_rate = monod_sing(C, C0, Rmax, Km, Km0)
%MONOD Summary of this function goes here
%   Detailed explanation goes here
    
    growth_rate = [0 0]';
    
    growth_rate(1) = (Rmax(1).*(C0./ (Km0(1) + C0)));
    growth_rate(2) = (Rmax(2).*C./ (Km + C)) .* (C0./ (Km0(2) + C0));
end

