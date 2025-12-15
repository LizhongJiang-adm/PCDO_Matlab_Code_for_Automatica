function [g_max_per_interval, t_max_per_interval] = MaxPathCstr_SubIntvlV3(decisionVars, constraintGrid, problem)
% =========================================================================
% Find Path Constraint Maximums Over Specified Sub-Intervals (V3)
% =========================================================================
% Description:
% This function simulates the system trajectory and finds the maximum value
% of a specific path constraint within each sub-interval defined by the
% `constraintGrid`. It is a core component for checking algorithm
% convergence and refining the constraint grid.
% This version is self-contained and does not rely on global variables.
%
% Inputs:
%   decisionVars   - The vector of decision variables (u).
%   constraintGrid - A vector of time points defining the sub-intervals
%                    for which to find the maximum constraint value.
%   problem        - The problem definition struct.
%
% Outputs:
%   g_max_per_interval - A vector containing the max value of the path
%                        constraint in each sub-interval.
%   t_max_per_interval - A vector containing the time at which each max
%                        occurs.
% =========================================================================

%% --- 1. Extract Parameters from Problem Struct ---
x0 = problem.state.x0;
controlGrid = problem.control.grid;
numDecisionVars = problem.numVars;
odeFun = problem.functions.odeSensitivity;
% Note: Assumes we are always checking the first path constraint
pathConstraintHandle = problem.pathConstraints.handle{1};

numStates = numel(x0);
numIntervals = problem.control.N;
numSensitivityStates = numStates * numDecisionVars;

%% --- 2. Reshape Decision Variables ---
controlParams = reshape(decisionVars, 1, numIntervals);

%% --- 3. Simulate Trajectory Over a Combined Grid ---
% The simulation needs to stop at every point in both the control grid
% and the constraint grid to ensure accurate evaluation.
combinedGrid = unique([controlGrid, constraintGrid]);
numCombinedIntervals = numel(combinedGrid) - 1;

numConstraintIntervals = numel(constraintGrid) - 1;
g_max_per_interval = zeros(1, numConstraintIntervals);
t_max_per_interval = zeros(1, numConstraintIntervals);

% Initialize the full state vector [states, sensitivities]
y0 = [x0, zeros(1, numSensitivityStates)];

% Store the full time history to easily find values in sub-intervals
full_time_history = [];
full_state_history = [];
full_control_history = [];

opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

for i = 1:numCombinedIntervals
    t_start = combinedGrid(i);
    t_end = combinedGrid(i+1);
    
    % Find which control interval this simulation segment belongs to
    % The `find` will get the index of the control interval that starts at or before t_start
    currentControlIdx = find(controlGrid <= t_start + 1e-10, 1, 'last');
    currentControlParams = controlParams(:, currentControlIdx);
    
    % Define the ODE for this small segment
    ode_interval = @(t, y) odeFun(t, y, currentControlParams, currentControlIdx, numIntervals, controlGrid(currentControlIdx:currentControlIdx+1));
    
    [t, y_full] = ode45(ode_interval, [t_start, t_end], y0);
    y0 = y_full(end, :); % Update initial condition for the next segment

    % Append results, avoiding duplicates
    if isempty(full_time_history)
        full_time_history = t;
        full_state_history = y_full;
    else
        full_time_history = [full_time_history; t(2:end)];
        full_state_history = [full_state_history; y_full(2:end,:)];
    end
end

% Extract only the state part from the full history
state_history = full_state_history(:, 1:numStates);

% Reconstruct the full control history
for i = 1:length(full_time_history)
    t_current = full_time_history(i);
    currentControlIdx = find(controlGrid <= t_current + 1e-10, 1, 'last');
    % Clamp index to be safe at the very end
    if t_current == controlGrid(end)
        currentControlIdx = numIntervals;
    end
    currentControlParams = controlParams(:, currentControlIdx);
    intervalDuration = controlGrid(currentControlIdx+1) - controlGrid(currentControlIdx);
    relativeTime = t_current - controlGrid(currentControlIdx);
    % full_control_history(i,1) = GetCtrlValV2(currentControlParams, relativeTime, intervalDuration, controlClass);
    full_control_history(i,1) = currentControlParams*ones(numel(relativeTime),1);
end

%% --- 4. Find Maximums in Each Constraint Sub-Interval ---
for i = 1:numConstraintIntervals
    t_interval_start = constraintGrid(i);
    t_interval_end = constraintGrid(i+1);
    
    % Find the indices corresponding to the current constraint interval
    idx_in_interval = find(full_time_history >= t_interval_start & full_time_history <= t_interval_end);
    
    if isempty(idx_in_interval)
        % This can happen if an interval is extremely small.
        % In this case, we can evaluate the constraint at the start point.
        g_max_per_interval(i) = pathConstraintHandle(t_interval_start, state_history(full_time_history == t_interval_start, :), full_control_history(full_time_history == t_interval_start, :));
        t_max_per_interval(i) = t_interval_start;
        continue;
    end
    
    % Extract the relevant portions of the trajectory
    t_subset = full_time_history(idx_in_interval);
    state_subset = state_history(idx_in_interval, :);
    control_subset = full_control_history(idx_in_interval, :);
    
    % Evaluate the path constraint over this subset
    g_values_in_interval = pathConstraintHandle(t_subset, state_subset, control_subset);
    
    % Find the maximum value and its corresponding time
    [max_g, idx_max] = max(g_values_in_interval);
    
    g_max_per_interval(i) = max_g;
    t_max_per_interval(i) = t_subset(idx_max);
end

end