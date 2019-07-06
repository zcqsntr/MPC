

true_params = {1, 0.5, [480000., 480000.], [520000., 520000.], [0.6, 0.6], [0.00048776, 0.00000102115], [0.00006845928, 0.00006845928], [-0.0001, -0.0001; -0.0001, -0.0001]};

target_N = [250, 550]';
initial_C = find_initial_guess(target_N, true_params);

disp(initial_C)

x = [target_N; initial_C(3:5)];
u = initial_C(1:2);

disp(chemostat_derivatives_doub(x, u, true_params))
