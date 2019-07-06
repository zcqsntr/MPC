function growth_rate = monod_doub(C, C0, Rmax, Km, Km0)
% matches python implementation
    growth_rate = ((Rmax.*C)./ (Km + C)) .* (C0./ (Km0 + C0));
  
end

