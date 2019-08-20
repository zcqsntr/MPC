fig = openfig('one_min_sampling_no_ramp_new_target.fig');

axObjs = fig.Children;
dataObjs = axObjs.Children;
disp(axObjs)
disp(dataObjs)


N2 = axObjs(10).Children.YData;
N1 = axObjs(12).Children.YData;

% i = find(N1 < 6, 1, 'first');
% N2 = N2(1:i);
% N1 = N1(1:i);

figure
plot(N1, 'linewidth', 2)
hold on
plot(N2, 'linewidth', 2)
xlabel('Time (minutes)')
ylabel('Population (10^6 cells L^{-1})')
legend("N1", "N2")
