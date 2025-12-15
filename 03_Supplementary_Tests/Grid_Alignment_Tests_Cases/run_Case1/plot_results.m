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
[g_max, t_max, state_and_g, time_vector, ~] = PathCnstrProfileV3(solution.decisionVariables, problem);

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

for j = 1:numel(solution.finalConstraintGrid)
    xline(h_main_axes, solution.finalConstraintGrid(j), '--', 'Color', [0.7, 0.7, 0.7], ...
          'HandleVisibility', 'off', 'LineWidth', 0.5);
end

% Adjust y-limits
y_max_val = max(g_max);
tmp = ylim(h_main_axes);
y_min_val = tmp(1);
range = y_max_val - y_min_val;
ylim(h_main_axes, [y_min_val - 0.1*range, y_max_val + 0.5*range]);

% Finalize main plot aesthetics
xlabel(h_main_axes, '$t$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel(h_main_axes, '$g(u,t)$', 'Interpreter', 'latex', 'FontSize', 15);
% title(h_main_axes, 'Path Constraint Profile');
box(h_main_axes, 'on');
hold(h_main_axes, 'off');

%% --- 3. Create Subplot (Zoomed-in View) ---
SubPicPst = [0.281047621150624,0.238095238095242,0.295738093135089,0.169270320681385];
h_sub_axes = axes('Position', SubPicPst);
hold(h_sub_axes, 'on');

for i = 1:numPathConstraints
    plot(h_sub_axes, time_vector, g_values(:, i), 'LineWidth', 2, 'HandleVisibility', 'off');
    if g_max(i) > 0
        plot(h_sub_axes, t_max(i), g_max(i), 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'HandleVisibility', 'off');
    else
        plot(h_sub_axes, t_max(i), g_max(i), 'bp', 'MarkerFaceColor', 'b', 'MarkerSize', 8, 'HandleVisibility', 'off');
    end
end
plot(h_sub_axes, time_vector, zeros(size(time_vector)), 'm--', 'HandleVisibility', 'off');

% Adjust limits for zoom
y_max_zoom = 0; y_min_zoom = g_max(1); % Assuming we zoom on the first constraint's max
range_zoom = y_max_zoom - y_min_zoom;
ylim(h_sub_axes, [y_min_zoom - 0.35*range_zoom, y_max_zoom + 0.35*range_zoom]);
xlim(h_sub_axes, [t_max(1) - 1e-5, t_max(1) + 1e-5]);
box(h_sub_axes, 'on');
hold(h_sub_axes, 'off');

%% --- 4. Draw Arrow from Data Point to Subplot ---
% =========================================================================
%                        ARROW DRAWING LOGIC
% =========================================================================

% --- Step A: Get main axes position and data limits ---
main_axes_pos = get(h_main_axes, 'Position');
main_axes_xlim = get(h_main_axes, 'XLim');
main_axes_ylim = get(h_main_axes, 'YLim');

% --- Step B: Convert the data point [t_max, g_max] to Figure coordinates ---
% We assume we are pointing from the first (or only) max point
x_point = t_max(1);
y_point = g_max(1);

norm_x_in_axes = (x_point - main_axes_xlim(1)) / (main_axes_xlim(2) - main_axes_xlim(1));
norm_y_in_axes = (y_point - main_axes_ylim(1)) / (main_axes_ylim(2) - main_axes_ylim(1));

point_x_fig = main_axes_pos(1) + norm_x_in_axes * main_axes_pos(3);
point_y_fig = main_axes_pos(2) + norm_y_in_axes * main_axes_pos(4);

% --- Step C: Define arrow start and end points in Figure coordinates ---
arrow_start = [point_x_fig, point_y_fig];

% Arrow end point: center of the subplot
arrow_end = [SubPicPst(1) + SubPicPst(3)/2, SubPicPst(2) + SubPicPst(4)/2];

arrow_vec =  arrow_end - arrow_start;
arrow_end = arrow_end - 0.3*arrow_vec;
arrow_start = arrow_start + 0.2*arrow_vec;
% --- Step D: Draw the annotation ---
annotation('textarrow', [arrow_start(1), arrow_end(1)], ...
                       [arrow_start(2), arrow_end(2)], ...
                        'FontSize', 10);
                       
% % You can also draw a box around the zoom area on the main plot
% zoom_box_x = [t_max(1)-1e-5, t_max(1)+1e-5, t_max(1)+1e-5, t_max(1)-1e-5, t_max(1)-1e-5];
% zoom_box_y_range = ylim(h_sub_axes);
% zoom_box_y = [zoom_box_y_range(1), zoom_box_y_range(1), zoom_box_y_range(2), zoom_box_y_range(2), zoom_box_y_range(1)];
% plot(h_main_axes, zoom_box_x, zoom_box_y, 'r--'); % Draw a red dashed box

%% --- 5. Finalize Legend ---
% Place the legend at the bottom of the figure
legend(h_main_axes, 'show', 'Interpreter', 'latex', 'FontSize', 15);

fprintf('Plotting complete.\n');
end