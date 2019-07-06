function dx = chemostat_derivatives_doub(x, u, params, is_measure)
% all derivatives exaclty the same as python implementation for the first
% step of double auxotroph
  
    N = x(1:2);
    C = x(3:4);
    C0 = x(5);
    
    C0in = params(1); q = params(2); y = params(3); y3 = params(4); Rmax = params(5); Km = params(6); Km0 = params(7); 
    C0in = cell2mat(C0in); q = cell2mat(q); y=cell2mat(y).'; y3=cell2mat(y3).'; Rmax = cell2mat(Rmax).'; Km = cell2mat(Km).'; Km0 = cell2mat(Km0).';
    
   
    growth_rates = monod_doub(C, C0, Rmax, Km, Km0);

    dN = N.*(growth_rates - q);
      
 
    dC = q.*(u - C) - (1./y).*growth_rates.*N;
  
   
    %disp(q.*(u - C) ); disp(- (1./y).*growth_rates.*N);
    %disp(u);disp(size(u));disp(size(C));disp(size(y));disp(size(growth_rates));disp(size(N));
       
    dC0 = q.*(C0in - C0) - sum(1./y3.*growth_rates.*N);
    
    
    % see if this is a measurement of the system or a simulation of the state model
    if is_measure == 1
        dx = dN;
    else 
        
        dx = [dN;dC;dC0];
    end
    
end

