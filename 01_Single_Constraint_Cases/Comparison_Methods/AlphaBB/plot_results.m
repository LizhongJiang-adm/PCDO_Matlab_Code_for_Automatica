function plot_results(solution, problem, stats)
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

% Recalculate the full path constraint profile using the refactored function
[g_max, t_max, state_and_g, time_vector, ~] = PathCnstrProfileV3(solution.decisionVariables, problem);

%% --- 2. Replicate the Original Plotting Logic ---
figure('Name', 'Path Constraint Profile (Publication Style)', 'NumberTitle', 'off');
hold on;

numStates = numel(problem.state.x0);
g_values = state_and_g(:, numStates+1:end); % Extract constraint values

numPathConstraints = problem.pathConstraints.num;

for i = 1:numPathConstraints
    % Plot the path constraint profile line
    plot(time_vector, g_values(:, i), 'LineWidth', 1.2, 'DisplayName', sprintf('$g_{%d}(t)$', i));

    % Plot the maximum point marker with conditional coloring and legend
    if g_max(i) > 0
        lgndGmax = ['$\max\limits_{t\in T} g_{', num2str(i), '}(t)>0$'];
        plot(t_max(i), g_max(i), 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', lgndGmax);
    else
        lgndGmax = ['$\max\limits_{t\in T} g_{', num2str(i), '}(t) \le 0$'];
        plot(t_max(i), g_max(i), 'bp', 'MarkerFaceColor', 'b', 'MarkerSize', 8, 'DisplayName', lgndGmax);
    end
end

% Plot the zero line (g=0 boundary)
plot(time_vector, zeros(size(time_vector)), 'm--', 'DisplayName', '$g(t)=0$');

% Plot the final constraint grid points from the solver
% for j = 1:numel(solution.finalConstraintGrid)
%     xline(solution.finalConstraintGrid(j), '--', 'Color', [0.7, 0.7, 0.7], ...
%           'HandleVisibility', 'off', 'LineWidth', 0.5);
% end

hold off;

%% --- 3. Finalize Plot Aesthetics ---
% Set labels and title with LaTeX interpreter
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$g(t)$', 'Interpreter', 'latex');
% title(sprintf('Path Constraint Profile\n(Iterations: %d, CPU Time: %.2f s)', ...
%       stats.totalIterations, stats.cpuTime));
title(sprintf('Path Constraint Profile'));


% Set legend with LaTeX interpreter
legend('show', 'Interpreter', 'latex', 'FontSize', 11, 'Location', 'northeast');

grid on;
box on;

fprintf('Plotting complete.\n');

end