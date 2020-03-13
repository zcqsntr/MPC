
clear;
params = {1, 0.5, [480000., 480000.], [520000., 520000.], [1., 1.1], [0.00048776, 0.00000102115], [0.00006845928, 0.00006845928], [0, 0; 0, 0]};

duration = 100;


options = odeset('NonNegative', [1 2 3 4 5]);


%% RUN WHOLE TIMESERIES AT ONCE, THIS DIVERGES

u_on = 0.1; % the on state of the controller
%odefun = @(t, x) chemostat_derivatives_doub(x, u, params, 0);
%[t_out, x_out] = ode45(odefun, [0 duration], x, options);

%% RUN STEP BY STEP, THIS OSCILLATES 

all_n_mins = [1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60];

for j = 1:length(all_n_mins)
    n_mins = all_n_mins(j);
    disp(n_mins);
    n_episodes = 30;
    Ts = n_mins/60; 
    tmax = round((24*60)/n_mins);
    trajectorys = {tmax,1};
    input_outputs = {tmax,1};
    in_outs = [];
    inputs = [];
    outputs = [];

    for i = 1:n_episodes
        disp(i);
        x0 = [20000; 30000; 0.; 0; 1];
        x = x0;
        xs = [x];
        uHistory = [];



        for t = 1:tmax
            u_in = randi([0 1], 2,1)*u_on;
            
            % add 10% pump noise 
            u = normrnd(u_in, u_in*0.1);
            
            uHistory = [uHistory u];
            odefun = @(t, x) chemostat_derivatives_doub(x, u, params);
            [t_out, x_out] = ode45(odefun, [0 Ts], x, options);

            x = x_out(end, :)';
            xs = [xs x];
            input = u_in;
            
            output = normrnd(x(1:2), x(1:2)*0.05); % add 5% measurement noise

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

    save("inputs_N1_"+n_mins+"_mins.mat", 'inputs_N1');
    save("inputs_N2_"+n_mins+"_mins.mat"', 'inputs_N2');
    save("outputs_N1_"+n_mins+"_mins.mat", 'outputs_N1');
    save("outputs_N2_"+n_mins+"_mins.mat", 'outputs_N2');
end
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