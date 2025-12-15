function [g_max_per_interval, t_max_per_interval] = MaxPathCstr_SubIntvlV4(decisionVars, constraintGrid, problem)
% =========================================================================
% Find Path Constraint Maximums Over Specified Sub-Intervals (V4 - Segmented Integration)
% =========================================================================
% Description:
% This function finds the maximum value of a path constraint within each
% sub-interval defined by 'constraintGrid'. This version uses a segmented
% integration approach, integrating the dynamics separately over each
% constraint sub-interval.
%
% ... (Inputs, Outputs are the same) ...
% =========================================================================

%% --- 1. Extract Parameters from Problem Struct ---
x0 = problem.state.x0;
controlGrid = problem.control.grid;
numDecisionVars = problem.numVars; % N
odeFun = problem.functions.odeSensitivity;
pathConstraintHandle = problem.pathConstraints.handle{1};

numStates = numel(x0);
numIntervals = problem.control.N;
numSensitivityStates = numStates * numDecisionVars;

%% --- 2. Initialize ---
numConstraintIntervals = numel(constraintGrid) - 1;
g_max_per_interval = zeros(1, numConstraintIntervals);
t_max_per_interval = zeros(1, numConstraintIntervals);

% Initialize the full state vector for integration
y0 = [x0, zeros(1, numSensitivityStates)];
opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

% Reshape decision vars for easier access
controlParams = reshape(decisionVars, 1, numIntervals);

%% --- 3. Segmented Integration and Analysis Loop ---
for i = 1:numConstraintIntervals
    t_start = constraintGrid(i);
    t_end = constraintGrid(i+1);
    
    % --- A. Integrate the trajectory ONLY over the current constraint interval ---
    
    % Find which control is active at the beginning of this interval
    currentControlIdx = find(controlGrid <= t_start + 1e-10, 1, 'last');
    currentControl = controlParams(currentControlIdx);
    
    % Check if the control switches INSIDE this constraint interval
    control_switch_point = controlGrid(currentControlIdx + 1);
    
    if control_switch_point > t_start && control_switch_point < t_end
        % --- Sub-case: Control switches within the interval ---
        % Integrate in two parts
        
        % Part 1: [t_start, control_switch_point]
        ode_h1 = @(t,y) odeFun(t, y, currentControl, currentControlIdx, numIntervals, []);
        [t1, y1] = ode45(ode_h1, [t_start, control_switch_point], y0, opts);
        
        % Part 2: [control_switch_point, t_end]
        nextControlIdx = currentControlIdx + 1;
        nextControl = controlParams(nextControlIdx);
        ode_h2 = @(t,y) odeFun(t, y, nextControl, nextControlIdx, numIntervals, []);
        [t2, y2] = ode45(ode_h2, [control_switch_point, t_end], y1(end,:), opts);
        
        % Concatenate results for this interval
        tout_interval = [t1; t2(2:end)];
        xout_interval = [y1; y2(2:end, :)];
        
    else
        % --- Sub-case: No control switch within the interval ---
        % Integrate in one part
        ode_handle = @(t,y) odeFun(t, y, currentControl, currentControlIdx, numIntervals, []);
        [tout_interval, xout_interval] = ode45(ode_handle, [t_start, t_end], y0, opts);
    end
    
    % --- B. Find the maximum of the path constraint in this interval ---
    tout_interval = tout_interval(2:end-1);
    state_subset = xout_interval(2:end-1, 1:numStates);
    
    % Reconstruct control for this interval's trajectory
    control_subset = zeros(size(tout_interval));
    for k = 1:length(tout_interval)
        ctrl_idx = find(controlGrid <= tout_interval(k), 1, 'last');
        if tout_interval(k) == controlGrid(end), ctrl_idx = numIntervals; end
        control_subset(k) = controlParams(ctrl_idx);
    end
    
    g_values_in_interval = pathConstraintHandle(tout_interval, state_subset, control_subset);
    [max_g, idx_max] = max(g_values_in_interval);
    
    g_max_per_interval(i) = max_g;
    t_max_current = tout_interval(idx_max);
    
    % --- C. Apply the nudge logic ---
    nudge = 1e-5;
    if abs(t_max_current - t_start) < 1e-9
        t_max_current = t_max_current + nudge;
    elseif abs(t_max_current - t_end) < 1e-9
        t_max_current = t_max_current - nudge;
    end
    
    t_max_per_interval(i) = max(t_start, min(t_end, t_max_current));
    
    % --- D. Update initial state for the NEXT constraint interval ---
    y0 = xout_interval(end, :);
end

end