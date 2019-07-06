% two species chemostat has 5 states
clear;
nx = 5; % number of states, should this be 7 to include the Cins 
ny = 2; % number of outputs 
nu = 2; % number of inputs

nl_mpc = nlmpc(nx, ny, nu);

Ts = 1; % sample time for prediction model
nl_mpc.Ts = Ts; % sample time for controller

nl_mpc.PredictionHorizon = 10;
nl_mpc.ControlHorizon = 4;


real_params = {0.5, 0.5, [480000000000., 480000000000.], [520000000000., 520000000000.], [0.6, 0.7], [0.00049, 0.00000102115], [0.00006845928, 0.00006845928]};
estimated_params = {0.5, 0.5, [480000000000., 480000000000.], [520000000000., 520000000000.], [0.6, 0.7], [0.00049, 0.00000102115], [0.00006845928, 0.00006845928]};


nl_mpc.Model.StateFcn = @(x,u) chemostat_derivatives_doub(x, u, estimated_params, 0);
nl_mpc.Model.OutputFcn = @(x,u) chemostat_derivatives_doub(x, u, real_params, 1);
% to improve efficiency specify jacobians of these functions

nl_mpc.MV(1).Min = 0;
nl_mpc.MV(1).Max = 0.3;

nl_mpc.MV(2).Min = 0;
nl_mpc.MV(2).Max = 0.3;

nl_mpc.Model.NumberOfParameters = 0;
% for the cost function (find out what this is exactly) or make a custm
% cost function
nl_mpc.Weights.OutputVariables = [0.000000001,0.000000001]; % penalisez distance from desired output
nl_mpc.Weights.ManipulatedVariablesRate = [0, 0]; % penalises large changes in control actions
nl_mpc.Weights.ManipulatedVariables = [0, 0]; % penalises distance from disired input
%nl_mpc.Optimization.CustomCostFcn = @(X, U, data, i) cost(X(1:2));
%nl_mpc.Optimization.ReplaceStandardCost = true;
x0 = [50000000000, 50000000000, 0.1, 0.1, 0.5]';

u0 = [0.15, 0.15]';


validateFcns(nl_mpc, x0, u0, []);


% set max iters for the MPC

%nl_mpc.Optimization.SolverOptions.MaxIter = 50;
%% Closed loop simulation


%state estimation for the hidden states
StateFcn = @(xk, uk) xk + chemostat_derivatives_doub(xk, uk, estimated_params, 0) * Ts; % UKF requires the state of system at time k+1, use FE for now
MeasureFcn = @(xk) xk(1:2);
UKF = unscentedKalmanFilter(StateFcn, MeasureFcn, x0); % use unscented over extended kalman filter

x_sys = x0;

UKF.State = x_sys;

uk = u0;

nloptions = nlmpcmoveopt;
nloptions.Parameters = {Ts};

Duration = 10;
hbar = waitbar(0, 'simulation progress');
xHistory = x_sys;
costs = [];
us = [];
y_targ = zeros(1, nl_mpc.PredictionHorizon);
y_targ(:,1) = 50000000000;
y_targ(:,2) = 50000000000;

% solver options
opts = odeset('NonNegative', [1 2 3 4 5]);

disp('ENTERING LOOP')
for ct = 1:(Duration/Ts)
   y = x_sys(1:2); %+ randn(2,1)*0.01;
   % use kalman filter to improve estimate 
   x_est = x_sys;%correct(UKF, y); % correct with kalman filter
   %compute optimal control moves
   
   [uk, nloptions, info] = nlmpcmove(nl_mpc, x_est, uk, y_targ);
   costs = [costs info.Cost];

   %predict  model states for next iter 
   predict(UKF, uk);
   us = [us uk];
   disp(uk)
   % IMPLEMENT FIRST OPTIMAL CONTROL MOVE AND UPDATE PLANT STATES
   odefun = @(t, xk) chemostat_derivatives_doub(xk, uk, real_params, 0);
   
   [t_out, x_out] = ode45(odefun, [0 Ts], x_sys', opts);
   x_sys = x_out(end,:);
   % generate sensor data with some noise 
  
   xHistory = [xHistory x_out(end,:)'];
   waitbar(ct*Ts/Duration, hbar);
end

close(hbar);

%% plot results

figure
subplot(3,2,1)
disp(size(xHistory(1,:)))
disp(size(0:Ts:Duration))
l1 = plot(xHistory(1,2:end));
legend(["N1"]);
xlabel('time')
ylabel('pop')
title('populations')

subplot(3,2,2)
disp(size(xHistory(2,:)))
disp(size(0:Ts:Duration))
l1 = plot(xHistory(2,2:end));
legend(["N2"]);
xlabel('time')
ylabel('pop')
title('populations')

subplot(3,2,3)
l1 = plot(xHistory(3,2:end));
hold on;
l2 = plot(xHistory(4,:));
legend(["AA1", "AA2"])
xlabel('time')
ylabel('conc')
title('Amino acids')


subplot(3,2,4)
plot(xHistory(5,2:end))
xlabel('time')
ylabel('conc')
title('Carbon')

subplot(3,2,5)
plot(costs)
xlabel('time')
ylabel('cost')
title('cost')

subplot(3,2,6)
plot(us')
xlabel('time')
ylabel('action')
title('actions')