function plot_constraint_composition(solution, problem)
% =========================================================================
% Plot the Composition of a Path Constraint
% =========================================================================
% Description:
% This function visualizes how a path constraint g(t) = u(t) + x1(t)/6
% is composed of its constituent parts: the control signal u(t) and the
% scaled state trajectory x1(t)/6.
%
% Inputs:
%   solution - The final solution struct from the solver.
%   problem  - The problem definition struct.
% =========================================================================

%% --- 1. Recalculate Full Trajectory ---
fprintf('Generating trajectories for constraint composition plot...\n');

% Get the full, dense trajectories for states, controls, and the path constraint
[g_max, ~, state_and_g, time_vector, control_vector] = ...
    PathCnstrProfileV3(solution.decisionVariables, problem);

% Extract the individual components
g_trajectory = state_and_g(:, end); % The full g(t) trajectory
x1_trajectory = state_and_g(:, 1);  % The x1(t) trajectory
u_trajectory = control_vector;      % The u(t) trajectory

%% --- 2. Calculate Constituent Parts of the Constraint ---
% Calculate the scaled state component of the constraint
x1_scaled_trajectory = x1_trajectory / 6;

%% --- 3. Create the Composition Plot ---
figure('Name', 'Path Constraint Composition', 'Position', [100, 100, 800, 500]);
hold on;

% --- Plot the three main curves ---
% Plot 1: The final path constraint g(t)
plot(time_vector, g_trajectory, 'b-', 'LineWidth', 2.5, ...
    'DisplayName', 'Path Constraint $g(t) = u_{traj}(t) + x_1(t)/6$');

% Plot 2: The control signal u(t)
plot(time_vector, u_trajectory, 'k-', 'LineWidth', 1.5, ...
    'DisplayName', 'Control Signal $u_{traj}(t)$');

% Plot 3: The scaled state x1(t)/6
plot(time_vector, x1_scaled_trajectory, 'r-', 'LineWidth', 2, ...
    'DisplayName', 'Scaled State $x_1(t)/6$');

% --- Add reference lines and markers ---
% Plot the g=0 boundary
yline(0, 'k--', 'HandleVisibility', 'off');

% Mark the maximum point of the path constraint
result_string = sprintf('$\\max_{t\\in[0,4.5]} g(t)=%.4e$', g_max);

[~, idx_max] = max(g_trajectory);
plot(time_vector(idx_max), g_trajectory(idx_max), 'p', ...
     'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', ...
     'MarkerSize', 10, 'DisplayName', result_string);

%% --- 4. Finalize Plot Aesthetics ---
hold off;
% title('Composition of the Path Constraint $g(t) = u(t) + x_1(t)/6$');
xlabel('Time (t)');
ylabel('Value');
legend('show', 'Interpreter', 'latex', 'Location', 'best',fontsize = 12);
grid on;
box on;

fprintf('Constraint composition plot generated successfully.\n');

end