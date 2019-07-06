
clear;
params = {0.5, 0.5, 480000000000., 520000000000., 0.6, 0.00049, 0.00006845928};
duration = 200;
x = [5e10, 0.1, 0.5]'; % initial X
xs = [x];

options = odeset('NonNegative', [1 2 3]);


%% RUN WHOLE TIMESERIES AT ONCE, THIS DIVERGES

u = [0.11]';
odefun = @(t, x) chemostat_derivatives_simp(x, u, params, 0);
[t_out, x_out] = ode45(odefun, [0 duration], x, options);

%% RUN STEP BY STEP, THIS OSCILLATES 
Ts = 1.; % reduce this to reduce amplitude of oscillations 
for t = 1:duration/Ts
    u = [0.3]'; % action
    disp(t)
    odefun = @(t, x) chemostat_derivatives_simp(x, u, params, 0);
    [t_out, x_out] = ode45(odefun, [0 Ts], x, options);
    
    x = x_out(end, :)';
 
    xs = [xs x];
end 

t_out = 1:length(xs)';
x_out = xs';


%% PLOT SOLUTIONS


%disp(chemostat_derivatives_doub(x0, u, params))
%disp("-----")
%disp(size(t_out))
%disp(size(x_ou

figure
subplot(3,2,1)
plot(t_out, x_out(:,1))
xlabel('time')
ylabel('pop')
title('N1 population')

subplot(3,2,2)
plot(t_out, x_out(:,2))
xlabel('time')
ylabel('conc')
title('C1')

subplot(3, 2,3)
plot(t_out, x_out(:,3))
xlabel('time')
ylabel('conc')
title('carbon')




%%
close all