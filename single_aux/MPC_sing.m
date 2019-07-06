
clc;
nx = 4; % number of states, should this be 7 to include the Cins 
ny = 2; % number of outputs 
nu = 1; % number of inputs

nlobj = nlmpc(nx,ny,nu);

params = {1., 0.4,[480000000000], [520000000000., 520000000000.], [0.6,0.7], 0.00000102115, [0.00006845928, 0.00006845928]};
nlobj.Model.StateFcn = @(x,u) chemostat_derivatives_sing(x, u , params);
nlobj.Model.OutputFcn = @(x,u) x(1:2);

nlobj_tracking = nlmpc(nx,ny,nu);

nlobj_tracking.Model.StateFcn = nlobj.Model.StateFcn;
nlobj_tracking.Model.OutputFcn = nlobj.Model.OutputFcn;
%nlobj_tracking.Jacobian.StateFcn = nlobj.Jacobian.StateFcn;

Ts = 5;
nlobj_tracking.Ts = Ts;
nlobj_tracking.PredictionHorizon = 10;
nlobj_tracking.ControlHorizon = 4;

nl_mpc.Weights.OutputVariables = [1]; % penalisez distance from desired output
nl_mpc.Weights.ManipulatedVariablesRate = [0]; % penalises large changes in control actions
nl_mpc.Weights.ManipulatedVariables = [0]; % penalises distance from desired input
nl_mpc.Weights.ECR = 0.000001; % slack variable tuning weight, larger this value the more likely controller is to violate constraint

for ct = 1:nu
    nlobj_tracking.MV(ct).Min = 0;
    nlobj_tracking.MV(ct).Max = 0.3;
end


%%
% Validate your prediction model and custom functions, and their Jacobians.
x0 = [50000000000, 50000000000, 0.1, 0.25]';
u0 = [0.11];
validateFcns(nlobj_tracking,x0,u0);

DStateFcn = @(xk,uk,Ts) chemostat_discrete_time(xk,uk,params, Ts);

DMeasFcn = @(xk) xk(1);

EKF = extendedKalmanFilter(DStateFcn,DMeasFcn,x0);
EKF.MeasurementNoise = 0.01;


Tsteps = 200;        
xHistory = x0';
uHistory = u0;
lastMV = u0;

Xref = zeros(nlobj_tracking.PredictionHorizon, 2);
Xref(:,1) = 10000000000;
Xref(:,2) = 20000000000;

% signal spaans for scale factors 
Uspan = 0.3;
Yspan = 100000000000;

nlobj_tracking.MV(1).ScaleFactor = Uspan;
nlobj_tracking.OV(1).ScaleFactor = Yspan;
opts = odeset('NonNegative', [1 2 3 4]);
hbar = waitbar(0,'Simulation Progress');
options = nlmpcmoveopt;
costs = [];
for k = 1:Tsteps
    % Obtain plant output measurements with sensor noise.
    yk = xHistory(k,1)' + randn*0.01;
    % Correct state estimation based on the measurements.
    xk = correct(EKF, yk);
    % Compute the control moves with reference previewing.
    [uk,options, info] = nlmpcmove(nlobj_tracking,xk,lastMV,Xref,[],options);
    costs = [costs info.Cost];
    % Predict the state for the next step.
    predict(EKF,uk,Ts);
    % Store the control move and update the last MV for the next step.
    uHistory(k,:) = uk';
    lastMV = uk;
    % Update the real plant states for the next step by solving the
    % continuous-time ODEs based on current states xk and input uk.
    ODEFUN = @(t,xk) chemostat_derivatives_sing(xk,uk, params);
    [TOUT,YOUT] = ode45(ODEFUN,[0 Ts], xHistory(k,:)', opts);
    % Store the state values.
    xHistory(k+1,:) = YOUT(end,:);            
    % Update the status bar.
    waitbar(k/Tsteps, hbar);
end
close(hbar)

%% 
% Compare the planned and actual closed-loop trajectories.
%% plot results

figure
subplot(3,2,1)

l1 = plot(xHistory(2:end, 1));
legend(["N1"]);
xlabel('time')
ylabel('pop')
title('population')

subplot(3,2,2)

l1 = plot(xHistory(2:end, 2));
legend(["N2"]);
xlabel('time')
ylabel('pop')
title('population')


subplot(3,2,3)
l1 = plot(xHistory(2:end, 3));
legend(["AA"])
xlabel('time')
ylabel('conc')
title('Amino acid')

subplot(3,2,4)
l1 = plot(xHistory(2:end, 4));
legend(["C"])
xlabel('time')
ylabel('conc')
title('Carbon')



subplot(3,2,5)
plot(costs)
xlabel('time')
ylabel('cost')
title('cost')

subplot(3,2,6)
plot(uHistory')
xlabel('time')
ylabel('action')
title('actions')

