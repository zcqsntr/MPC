function dx = chemostat_derivatives_sing(x, u, params)
%CHEMOSTAT_DERIVATIVES Summary of this function goes here
%   Detailed explanation goes here
    N = x(1:2);
    C = x(3);
    C0 = x(4);
    
    
    C0in= params(1); q= params(2); y= params(3); y0= params(4); Rmax= params(5); Km= params(6); Km0= params(7); 
    C0in = cell2mat(C0in); q = cell2mat(q); y=cell2mat(y).'; y0=cell2mat(y0).'; Rmax = cell2mat(Rmax).'; Km=cell2mat(Km).'; Km0 = cell2mat(Km0).';
    
   
    growth_rates = monod_sing(C, C0, Rmax, Km, Km0);
    
    
    dN = N.*(growth_rates - q);
    
    dC = q*(u - C) - (1/y)*growth_rates(2)*N(2);
    
    %disp("-----")
    %disp(growth_rates)
    %disp(dN)
    %disp(dC)
    %disp(u);disp(size(u));%disp(size(C));disp(size(y));disp(size(growth_rates));disp(size(N));
       
    dC0 = q.*(C0in - C0) - sum(1./y0.*growth_rates.*N);
    
    dx = [dN;dC;dC0];
  
    
end


