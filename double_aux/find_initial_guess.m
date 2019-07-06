


function [initial_C] = find_initial_guess(target_N, params)
%solves the nonlinear systems of ODEs to find an intial guess for Cin
    C0in = params(1); q = params(2); y = params(3); y3 = params(4); Rmax = params(5); Km = params(6); Km0 = params(7); A = params(8); 
    C0in = cell2mat(C0in); q = cell2mat(q); y=cell2mat(y).'; y3=cell2mat(y3).'; Rmax = cell2mat(Rmax).'; Km = cell2mat(Km).'; Km0 = cell2mat(Km0).'; A = cell2mat(A);
       
    N = target_N;
    
    x0 = [0.05,0.05, 0.1, 0.1, 0]';
    
    
    function F = chemostat_ders(all_Cs, target_N)
        %disp(q.*(u - C) ); disp(- (1./y).*growth_rates.*N);
        %disp(u);disp(size(u));disp(size(C));disp(size(y));disp(size(growth_rates));disp(size(N));
        Cins = all_Cs(1:2);
        Cs = all_Cs(3:4);
        C0 = all_Cs(5);
        
        
   
        growth_rates = monod_doub(Cs, C0, Rmax, Km, Km0);
        
        
        F(1:2) = N.*(growth_rates + A*N - q);
        
        
        F(3:4) = q.*(Cins - Cs) - (1./y).*growth_rates.*N;
        F(5) = q.*(C0in - C0) - sum(1./y3.*growth_rates.*N);

    end

    initial_C = fsolve(@(x)chemostat_ders(x), x0);
end



