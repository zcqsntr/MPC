function growth_rate = monod_sing(C0, Rmax, Km0)
%MONOD Summary of this function goes here
%   Detailed explanation goes here
    
    growth_rate = Rmax*(C0./ (Km0 + C0));
end

