function [g_profile_max, t_at_max, state_and_g, time_vector, control_vector] = PathCnstrProfileV3(decisionVars, problem)
% =========================================================================
% Calculate State Trajectory and Path Constraint Profile (Simplified Version)
% =========================================================================
% Description:
% This function simulates the system dynamics given a set of decision
% variables (control parameters). It computes the state trajectory, the
% control profile, and the value of the path constraint(s) over time.
% This version integrates the full state-sensitivity system for simplicity
% and then discards the sensitivity results.
%
% Inputs:
%   decisionVars - The vector of optimal decision variables (u).
%   problem      - The problem definition struct, containing all necessary
%                  parameters like initial state, dynamics, grid, etc.
%
% Outputs:
%   (Outputs are the same as before)
% =========================================================================

%% --- 1. Extract Parameters from the Problem Struct ---
x0 = problem.state.x0;
controlGrid = problem.control.grid;
numDecisionVars = problem.numVars;
odeFun = problem.functions.odeSensitivity;
pathConstraintHandles = problem.pathConstraints.handle;
numPathConstraints = problem.pathConstraints.num;

numStates = numel(x0);
numIntervals = problem.control.N;
numSensitivityStates = numStates * numDecisionVars;

%% --- 2. Reshape Decision Variables Based on Control Type ---

controlParams = reshape(decisionVars, 1, numIntervals);

%% --- 3. Simulate the System Trajectory ---
% %% REFAC: Your suggestion - integrate the full system for simplicity.
% Initialize the full state vector [states, sensitivities]
y0 = [x0, zeros(1, numSensitivityStates)]; 

% Initialize output arrays
time_vector = [];
full_state_out = []; % Store the full output from ODE solver
control_vector = [];

% Set high-precision ODE solver options
opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

% Integrate interval by interval
for i = 1:numIntervals
    t_start = controlGrid(i);
    t_end = controlGrid(i+1);
    time_span = [t_start, t_end];
    
    currentControlParams = controlParams(:, i);
    
    % The ODE function now directly matches what the solver expects. No wrapper needed.
    ode_interval = @(t, y) odeFun(t, y, currentControlParams, i, numIntervals, time_span);
    
    [t, y_full] = ode45(ode_interval, time_span, y0, opts);
    
    % Set the initial state for the next interval
    y0 = y_full(end, :);

    % Append results, avoiding duplicate time points between intervals
    if isempty(time_vector)
        time_vector = t;
        full_state_out = y_full;
    else
        time_vector = [time_vector; t(2:end)];
        full_state_out = [full_state_out; y_full(2:end, :)];
    end
    
    % Reconstruct the control history for plotting
    interval_duration = t_end - t_start;
    u_vals = currentControlParams*ones(numel(t - t_start),1);

    if isempty(control_vector)
        control_vector = u_vals;
    else
        control_vector = [control_vector; u_vals(2:end)];
    end
end

% %% REFAC: Extract only the state part (first numStates columns) for analysis.
state_out = full_state_out(:, 1:numStates);

%% --- 4. Evaluate Path Constraints Along the Trajectory ---
g_values = zeros(length(time_vector), numPathConstraints);
g_profile_max = zeros(1, numPathConstraints);
t_at_max = zeros(1, numPathConstraints);

for i = 1:numPathConstraints
    % Evaluate the i-th path constraint
    g_values(:, i) = pathConstraintHandles{i}(time_vector, state_out, control_vector);
    
    % Find its maximum value and the time it occurred
    [g_profile_max(i), idx_max] = max(g_values(:, i));
    t_at_max(i) = time_vector(idx_max);
end

%% --- 5. Package Final Outputs ---
state_and_g = [state_out, g_values];

end