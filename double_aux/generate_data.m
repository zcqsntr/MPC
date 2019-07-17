
clear;
params = {1, 0.5, [480000., 480000.], [520000., 520000.], [0.6, 0.6], [0.00048776, 0.00000102115], [0.00006845928, 0.00006845928], [-0.0001, -0.0001; -0.0001, -0.0001]};

duration = 100;


options = odeset('NonNegative', [1 2 3 4 5]);


%% RUN WHOLE TIMESERIES AT ONCE, THIS DIVERGES

u_on = 0.1; % the on state of the controller
%odefun = @(t, x) chemostat_derivatives_doub(x, u, params, 0);
%[t_out, x_out] = ode45(odefun, [0 duration], x, options);

%% RUN STEP BY STEP, THIS OSCILLATES 
Ts = 3; 
trajectorys = {50,1};
input_outputs = {50,1};
in_outs = [];
inputs = [];
outputs = [];

for i = 1:25
    x0 = [250; 550; 0.01472; 0.0013589; 0.009982];
    x = x0;
    xs = [x];
    uHistory = [];
    
    
    
    for t = 1:duration/Ts
        u = randi([0 1], 2,1)*u_on;
        uHistory = [uHistory u];
        odefun = @(t, x) chemostat_derivatives_doub(x, u, params);
        [t_out, x_out] = ode45(odefun, [0 Ts], x, options);

        x = x_out(end, :)';
        xs = [xs x];
        input = u;
        output = x(1:2);
        
        inputs = [inputs input];
        outputs = [outputs output];
        in_out = [input; output];
        in_outs = [in_outs in_out];
    
        if x(1) < 10 || x(2) < 10
            break
        end 

    end 
    trajectorys{i} = xs';
    
    % store input output 
    input_outputs{i} = in_outs';
end
inputs_N1 = inputs(1, :);
outputs_N1 = outputs(1,:);
inputs_N2 = inputs(2, :);
outputs_N2 = outputs(2,:);

t_out = 1:length(xs)';
Tsteps = t;
xHistory = xs';

%save('inputs_N1.mat', 'inputs_N1');
%save('inputs_N2.mat', 'inputs_N2');
%save('outputs_N1.mat', 'outputs_N1');
%save('outputs_N2.mat', 'outputs_N2');

%% PLOT SOLUTIONS


%disp(chemostat_derivatives_doub(x0, u, params))
%disp("-----")
%disp(size(t_out))
%disp(size(x_ou

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



subplot(4,2,7)
plot(uHistory')
xlabel('time')
ylabel('action')
title('actions')


%%
close all