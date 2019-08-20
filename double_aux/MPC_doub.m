clc;
clear;
nx = 5; % number of states, should this be 7 to include the Cins?
ny = 2; % number of outputs 
nu = 2; % number of inputs

nlobj = nlmpc(nx,ny,nu);

true_params = {1, 0.5, [480000., 480000.], [520000., 520000.], [0.6, 0.7], [0.00048776, 0.00000102115], [0.00006845928, 0.00006845928], [0.0, 0.0; 0.0, 0.0]};

est_params = {1, 0.5, [480000., 480000.], [520000., 520000.], [0.6, 0.7], [0.00048776, 0.00000102115], [0.00006845928, 0.00006845928], [0.0,0.0; 0.0, 0.0]};

nlobj.Model.StateFcn = @(x,u) chemostat_derivatives_doub(x, u , est_params);
nlobj.Model.OutputFcn = @(x,u) x(1:2); % add noise here if needed

%nlobj_tracking.Jacobian.StateFcn = nlobj.Jacobian.StateFcn;


Ts = 60/60;
nlobj.Ts = Ts;
nlobj.PredictionHorizon = 10;
nlobj.ControlHorizon = 1;

nl_mpc.Weights.OutputVariables = [1, 1]; % penalises distance from desired output
nl_mpc.Weights.ManipulatedVariablesRate = [0, 0]; % penalises large changes in control actions
nl_mpc.Weights.ManipulatedVariables = [0, 0]; % penalises distance from desired input
nl_mpc.Weights.ECR = 0.000001; % slack variable tuning weight, larger this value the more likely controller is to violate constraint

for ct = 1:nu
    nlobj.MV(ct).Min = 0;
    nlobj.MV(ct).Max = 0.1;
end


%%
% Validate your prediction model and custom functions, and their Jacobians.
target_N = [15000, 30000]';

x0 = [15000, 30000, 2.5, 0., 0.91]';
%x0 = [250, 550, 0.0, 0.0, 1.]'; % iniitial condition from the BioRXiv paper
xk = x0;
u0 = find_initial_guess(target_N, true_params);
u0 = u0(1:2)';
u0 = [0.034, 0.059];
disp(u0)

validateFcns(nlobj,x0,u0);

DStateFcn = @(xk,uk,Ts) chemostat_discrete_time(xk, uk, est_params, Ts);

DMeasFcn = @(xk) xk(1:2);

EKF = extendedKalmanFilter(DStateFcn,DMeasFcn,x0);
EKF.MeasurementNoise = 0.0;


Tsteps = 10000;        
xHistory = x0';
uHistory = u0;

lastMV = u0;


Xref = zeros(nlobj.PredictionHorizon, 2);
Xref(:,1) = target_N(1);
Xref(:,2) = target_N(2);

% signal spaans for scale factors 
Uspan = 0.1;
Yspan = 1000;

nlobj.MV(1).ScaleFactor = Uspan;
nlobj.MV(2).ScaleFactor = Uspan;

nlobj.OV(1).ScaleFactor = target_N(1);
nlobj.OV(2).ScaleFactor = target_N(2);

opts = odeset('NonNegative', [1 2 3 4 5]);
hbar = waitbar(0,'Simulation Progress');
options = nlmpcmoveopt;
costs = [];


t0 = 300;
t_col = 600;
increase = 0.1;


for k = 1:Tsteps
    % Obtain plant output measurements with sensor noise.
    yk = xHistory(k,1:2); % can add noise here
    % Correct state estimation based on the measurements.
    est_xk = correct(EKF, yk);
   % est_xk(est_xk<0) = 0; %kalman filter can somtimes output egatvie values 
    %est_xk = xHistory(k,:);
    % Compute the control moves with reference previewing.
   
    [uk,options, info] = nlmpcmove(nlobj,est_xk,lastMV,Xref,[],options);
    costs = [costs info.Cost];
    
    % Predict the state for the next step.
    predict(EKF,uk,Ts);
    % Store the control move and update the last MV for the next step.
    
    uHistory(k,:) = uk';
    lastMV = uk;
    % Update the real plant states for the next step by solving the
    % continuous-time ODEs based on current states xk and input uk.
    
%     if k > t0 && k < t0 + t_col
%         true_params{5}(1) = true_params{5}(1) + increase/t_col;
%     end
        
    ODEFUN = @(t,x) chemostat_derivatives_doub(x,uk, true_params);
    
    [TOUT,XOUT] = ode45(ODEFUN,[0 Ts], xHistory(k,:)', opts);
    % Store the state values.
    
    xHistory(k+1,:) = XOUT(end,:); 
    
    
   % xk = YOUT(end,:);
   
  
%     disp('-------------------')
%     disp(yk)
%     disp(xk(1:2))
    % Update the status bar.
    waitbar(k/Tsteps, hbar);
end
close(hbar)

%% 
% Compare the planned and actual closed-loop trajectories.
%% plot results

figure
subplot(4,2,1)
disp(size(xHistory(:,1)))
disp(size(0:Tsteps))
l1 = plot(xHistory(:, 1));
legend(["N1"]);
xlabel('time')
ylabel('pop')
title('population')

subplot(4,2,2)
disp(size(xHistory(:,2)))
disp(size(0:Tsteps))
l1 = plot(xHistory(:, 2));
legend(["N2"]);
xlabel('time')
ylabel('pop')
title('population')


subplot(4,2,3)
l1 = plot(xHistory(:, 3));
legend(["AA1"])
xlabel('time')
ylabel('conc')
title('Amino acid')

subplot(4,2,4)
l1 = plot(xHistory(:, 4));
legend(["AA2"])
xlabel('time')
ylabel('conc')
title('Amino acid')

subplot(4,2,5)
l1 = plot(xHistory(:, 5));
legend(["C0"])
xlabel('time')
ylabel('conc')
title('Carbon')



subplot(4,2,6)
plot(costs)
xlabel('time')
ylabel('cost')
title('cost')

subplot(4,2,7)
plot(uHistory)
xlabel('time')
ylabel('action')
title('actions')