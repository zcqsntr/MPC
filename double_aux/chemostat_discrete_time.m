function xk1 = chemostat_discrete_time(xk, uk, params, Ts)
    % Nonlinear discrete-time state transition for the flying robot model.
    %
    % Discretization uses the trapezoidal formula, assuming the input, uk, is
    % constant during time interval Ts.
    %
    % The trapezoidal formula is implicit; that is, you cannot solve for xk1
    % algebraically. You need to solve a set of nonlinear equations.  
    % This function uses FSOLVE for this purpose.

    % FlyingRobotStateFcn is the continuous-time state function for this
    % example. Obtain the state derivatives at the current point.
    
    ffun = @(xk,uk) chemostat_derivatives_doub(xk, uk, params);
    fk = ffun(xk,uk);
    
    % Extrapolation using xk1 = xk + Ts*fk is risky, since it might put xk1 in
    % an infeasible region, which could prevent convergence. A safer
    % alternative is xk1 = xk, but this method produces a poor estimate.
%     
%     xk1 = xk + Ts*fk;
%     
%     % Solve for xk1 satisfying the Trapezoidal rule.
%     FUN = @(xk1) TrapezoidalRule(xk,xk1,uk,Ts,fk,ffun);
%     Options = optimoptions('fsolve','Display','none');
%     xk1 = fsolve(FUN,xk1,Options);


    % solve using ode45]
    % xk, uk the right shape
    
   
    ODEFUN = @(t,x) chemostat_derivatives_doub(x,uk, params);
    %opts = odeset('NonNegative', [1 2 3 4 5]);
    [TOUT,XOUT] = ode45(ODEFUN,[0 Ts], xk);
    
    xk1 = XOUT(end,:)'; % (5,1)
   
    
end

% Trapezoidal rule function
function f = TrapezoidalRule(xk,xk1,uk,Ts,fk,ffun)
    % Time derivatives at point xk1.
    fk1 = ffun(xk1,uk);
    % The following must be zero to satisfy the Trapezoidal Rule
    f = xk1 - (xk + (Ts/2)*(fk1 + fk));
end

