function plot_pathconstraint(DecVar, problem)
% =========================================================================
% Plot Results for the Optimal Control Problem
% =========================================================================
% Description:
% This function visualizes the solution of the optimal control problem.
% It primarily plots the path constraint profile to show how the constraint
% is satisfied along the trajectory. It also shows the final refined grid
% points used by the upper-bound method.
%
% Inputs:
%   solution - A struct containing the final optimization results.
%   problem  - The problem definition struct.
%   stats    - A struct with statistics about the solution process.
% =========================================================================

%% --- 1. Extract Data and Recalculate Profile ---
fprintf('Generating path constraint profile for the final solution...\n');
[g_max, t_max, state_and_g, time_vector, ~] = PathCnstrProfileV3(DecVar, problem);

%% --- 2. Create Main Plot ---
figure('Name', 'Path Constraint Profile (Publication Style)', 'NumberTitle', 'off');

% --- Main Axes ---
% Store the handle to the main axes for coordinate conversion
h_main_axes = axes; % Adjusted position for legend below
hold(h_main_axes, 'on');

numStates = numel(problem.state.x0);
g_values = state_and_g(:, numStates+1:end);
numPathConstraints = problem.pathConstraints.num;

for i = 1:numPathConstraints
    plot(h_main_axes, time_vector, g_values(:, i), 'LineWidth', 2, 'DisplayName', sprintf('$g(u,t)$'));
    if g_max(i) > 0
        lgndGmax = ['$\max\limits_{t\in \Gamma} g(u,t)>0$'];
        plot(h_main_axes, t_max(i), g_max(i), 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', lgndGmax);
    else
        lgndGmax = ['$\max\limits_{t\in \Gamma} g(u,t) \le 0$'];
        plot(h_main_axes, t_max(i), g_max(i), 'bp', 'MarkerFaceColor', 'b', 'MarkerSize', 8, 'DisplayName', lgndGmax);
    end
end

plot(h_main_axes, time_vector, zeros(size(time_vector)), 'm--', 'DisplayName', '$g(u,t)=0$');

% for j = 1:numel(solution.finalConstraintGrid)
%     xline(h_main_axes, solution.finalConstraintGrid(j), '--', 'Color', [0.7, 0.7, 0.7], ...
%           'HandleVisibility', 'off', 'LineWidth', 0.5);
% end

% Adjust y-limits
y_max_val = 0;
tmp = ylim(h_main_axes);
y_min_val = tmp(1);
range = y_max_val - y_min_val;
ylim(h_main_axes, [y_min_val - 0.1*range, y_max_val + 0.4*range]);

% Finalize main plot aesthetics
xlabel(h_main_axes, '$t$', 'Interpreter', 'latex');
ylabel(h_main_axes, '$g(u,t)$', 'Interpreter', 'latex');
% title(h_main_axes, 'Path Constraint Profile');
box(h_main_axes, 'on');
hold(h_main_axes, 'off');

%% --- 5. Finalize Legend ---
% Place the legend at the bottom of the figure
legend(h_main_axes, 'show', 'Interpreter', 'latex', 'FontSize', 11);

fprintf('Plotting complete.\n');
end