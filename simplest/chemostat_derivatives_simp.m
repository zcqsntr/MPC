function dx = chemostat_derivatives_simp(x, u, params)
%CHEMOSTAT_DERIVATIVES Summary of this function goes here
%   Detailed explanation goes here
    N = x(1);
    C0 = x(2);
   
    
    q= params(1); y0= params(2); Rmax= params(3); Km0= params(4);
    q = cell2mat(q); y0=cell2mat(y0).'; Rmax = cell2mat(Rmax).';  Km0 = cell2mat(Km0).';
    
   
    growth_rate = monod_simp(C0, Rmax, Km0);
    
    
    dN = N.*(growth_rate - q);
  
    
    %disp("-----")
    %disp(growth_rates)
    %disp(dN)
    %disp(dC)
    %disp(u);disp(size(u));%disp(size(C));disp(size(y));disp(size(growth_rates));disp(size(N));
       
    dC0 = q.*(u - C0) - 1./y0.*growth_rate.*N;
    
    
    % see if this is a measurement of the system or a simulation of the state model
    
    dx = [dN;dC0];
    
end


