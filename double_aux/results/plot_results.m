fig = openfig('1000_tsteps_good_ICS.fig');

axObjs = fig.Children;
dataObjs = axObjs.Children;
disp(axObjs)
disp(dataObjs)


N2 = axObjs(10).Children.YData;
N1 = axObjs(12).Children.YData;



figure
plot(N1, 'linewidth', 2)
hold on
plot(N2, 'linewidth', 2)
xlabel('Time (minutes)')
ylabel('Population (10^6 cells L^{-1})')
legend("N1", "N2")
