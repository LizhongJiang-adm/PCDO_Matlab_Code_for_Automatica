% --- Post-processing and Visualization using Subplots ---
w_opt = solution.optimal_weights;
ocp = problem.ocp;
grid_points = problem.grid_points;

% Recreate strategy for helper methods
strategy = SequentialShootingStrategy_UB(ocp, grid_points);

% Get dense trajectory for plotting
[~, ~, ~, T_sol] = strategy.unpack_solution(w_opt);
t_dense = linspace(0, T_sol, 2000);
[x_dense, u_dense] = strategy.get_dense_trajectory(w_opt, t_dense);

% Create a CasADi function to evaluate all path constraints
path_fun = casadi.Function('PathFun', {ocp.model.x_sym, ocp.model.u_sym}, {vertcat(problem.path_constraints{:})});
path_values_dense = full(path_fun(x_dense, u_dense));

% --- Create Figure with Two Subplots ---
figure('Name', 'PFBF Path Constraint Trajectories', 'Position', [100, 100, 700, 600]);

% Define colors for consistency
color1 = [0, 0.4470, 0.7410]; % Blue
color2 = [0.8500, 0.3250, 0.0980]; % Orange-Red
magenta_dash = [0.9290, 0.6940, 0.1250]; % A standard magenta/pink color

% --- Subplot for Path Constraint 1 ---
ax1 = subplot(2, 1, 1);
plot(ax1, t_dense, path_values_dense(1,:), 'LineWidth', 1.5, 'Color', color1, 'DisplayName', 'Constraint Trajectory');
hold(ax1, 'on');
yline(ax1, 0, 'm--', 'DisplayName', 'Boundary g=0'); % MODIFICATION: Changed to 'm--'

% Find and plot the maximum value for this constraint
[g1_max, idx1_max] = max(path_values_dense(1,:));
t1_max = t_dense(idx1_max);
plot(ax1, t1_max, g1_max, 'p', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red', 'MarkerSize', 8, 'DisplayName', 'Maximum Value'); % MODIFICATION: Removed black edge

% MODIFICATION: Dynamic Y-axis limits
min_g1 = min(path_values_dense(1,:));
range1 = g1_max - min_g1;
ylim(ax1, [min_g1 - 0.2*range1, g1_max + 0.2*range1]);

% Dynamic Title with Max Value
title_str1 = sprintf('Path Constraint 1: $x_1 - 40 \\leq 0$ (Max Value: %+.4e)', g1_max);
title(ax1, title_str1, 'Interpreter', 'latex');
ylabel(ax1, 'Constraint Value $g_1(t)$', 'Interpreter', 'latex');
grid(ax1, 'on');
box(ax1, 'on');
legend(ax1, 'show', 'Location', 'best');
hold(ax1, 'off');

% --- Subplot for Path Constraint 2 ---
ax2 = subplot(2, 1, 2);
plot(ax2, t_dense, path_values_dense(2,:), 'LineWidth', 1.5, 'Color', color2, 'DisplayName', 'Constraint Trajectory');
hold(ax2, 'on');
yline(ax2, 0, 'm--', 'DisplayName', 'Boundary g=0'); % MODIFICATION: Changed to 'm--'

% Find and plot the maximum value for this constraint
[g2_max, idx2_max] = max(path_values_dense(2,:));
t2_max = t_dense(idx2_max);
plot(ax2, t2_max, g2_max, 'p', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red', 'MarkerSize', 8, 'DisplayName', 'Maximum Value'); % MODIFICATION: Removed black edge

% MODIFICATION: Dynamic Y-axis limits
min_g2 = min(path_values_dense(2,:));
range2 = g2_max - min_g2;
ylim(ax2, [min_g2 - 0.2*range2, g2_max + 0.2*range2]);

% Dynamic Title with Max Value
title_str2 = sprintf('Path Constraint 2: $x_2 - 0.5 \\leq 0$ (Max Value: %+.4e)', g2_max);
title(ax2, title_str2, 'Interpreter', 'latex');
xlabel(ax2, 'Time (s)');
ylabel(ax2, 'Constraint Value $g_2(t)$', 'Interpreter', 'latex');
grid(ax2, 'on');
box(ax2, 'on');
legend(ax2, 'show', 'Location', 'best');
hold(ax2, 'off');

% Link the X-axes so that zooming/panning in one affects the other
linkaxes([ax1, ax2], 'x');