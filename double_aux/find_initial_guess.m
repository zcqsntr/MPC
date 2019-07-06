function [initial_C] = find_initial_guess(target_N, params)
%solves the nonlinear systems of ODEs to find an intial guess for Cin
    C0in = params(1); q = params(2); y = params(3); y3 = params(4); Rmax = params(5); Km = params(6); Km0 = params(7); A = params(8); 
    C0in = cell2mat(C0in); q = cell2mat(q); y=cell2mat(y).'; y3=cell2mat(y3).'; Rmax = cell2mat(Rmax).'; Km = cell2mat(Km).'; Km0 = cell2mat(Km0).'; A = cell2mat(A);
       
    N = target_N;

    
    
    function F = chemostat_ders(Cs)


        %disp(q.*(u - C) ); disp(- (1./y).*growth_rates.*N);
        %disp(u);disp(size(u));disp(size(C));disp(size(y));disp(size(growth_rates));disp(size(N));
        
        growth_rates = monod_doub(Cs(1:2), Cs(3), Rmax, Km, Km0);
        F(1:2) = N.*(growth_rates + A*N - q);
        F(3:4) = q.*(u - Cs(1:2)) - (1./y).*growth_rates.*N;
        F(5) = q.*(C0in - Cs(3)) - sum(1./y3.*growth_rates.*N);
    
    end
    
    x0 = [0,0];
    initial_C = fsolve(chemostat_ders, x0);
    
    
    
    
    
   
end

