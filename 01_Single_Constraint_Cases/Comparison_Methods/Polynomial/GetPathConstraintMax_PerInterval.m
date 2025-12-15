function [g_max_per_interval, t_max_per_interval] = GetPathConstraintMax_PerInterval(decisionVars, problem, CstrIntvl)
% =========================================================================
% Find True Path Constraint Maximums Over Specified Sub-Intervals
% =========================================================================
% Description:
% This function simulates the system trajectory (including sensitivities) and
% finds the maximum value of the true path constraint within each sub-interval
% defined by 'CstrIntvl'. This version is self-contained and assumes constant control.
%
% Inputs:
%   decisionVars     - The vector of control variables (u).
%   problem          - The problem definition struct.
%   CstrIntvl        - A vector of time points defining the sub-intervals
%                      for which to find the maximum constraint value.
%
% Outputs:
%   g_max_per_interval - A vector containing the max value of the path
%                        constraint in each sub-interval.
%   t_max_per_interval - A vector containing the time at which each max occurs.
% =========================================================================

%% --- 1. Extract Parameters & Initialize ---
N = problem.control.N;
x0 = problem.state.x0;
controlGrid = problem.control.grid;
numStates = problem.system.dimensions;
odeFun = problem.functions.odeSensitivity;
pathConstraintHandle = problem.pathConstraints.handle{1};

num_cstr_intervals = length(CstrIntvl) - 1;
g_max_per_interval = zeros(1, num_cstr_intervals);
t_max_per_interval = zeros(1, num_cstr_intervals);

%% --- 2. Simulate Trajectory (with Sensitivities) Over a Combined Grid ---
combinedGrid = unique([controlGrid, CstrIntvl]);
num_combined_intervals = length(combinedGrid) - 1;

% Initialize the full state vector, including sensitivities
numSensitivityStates = numStates * N;
y0 = [x0, zeros(1, numSensitivityStates)];

full_time_history = [];
full_xout_history = []; % This will store the full state + sensitivity vector
% opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

for i = 1:num_combined_intervals
    t_start = combinedGrid(i);
    t_end = combinedGrid(i+1);
    
    control_idx = find(controlGrid <= t_start + 1e-10, 1, 'last');
    u_i = decisionVars(control_idx);
    
    ode_handle_interval = @(t, y) odeFun(t, y, u_i, control_idx, N, [t_start, t_end]);
    
    [t_interval, x_interval] = ode45(ode_handle_interval, [t_start, t_end], y0);
    y0 = x_interval(end, :);

    if isempty(full_time_history)
        full_time_history = t_interval;
        full_xout_history = x_interval;
    else
        full_time_history = [full_time_history; t_interval(2:end)];
        full_xout_history = [full_xout_history; x_interval(2:end,:)];
    end
end

%% --- 3. Reconstruct Control History & Evaluate Path Constraint ---
% Extract only the state part for constraint evaluation
full_state_history = full_xout_history(:, 1:numStates);

full_control_history = zeros(length(full_time_history), 1);
for i = 1:length(full_time_history)
    t_current = full_time_history(i);
    control_idx = find(controlGrid <= t_current + 1e-10, 1, 'last');
    if t_current == controlGrid(end), control_idx = N; end
    full_control_history(i) = decisionVars(control_idx);
end

full_g_values = pathConstraintHandle(full_time_history, full_state_history, full_control_history);

%% --- 4. Find Maximums in Each Constraint Sub-Interval ---
for i = 1:num_cstr_intervals
    t_interval_start = CstrIntvl(i);
    t_interval_end = CstrIntvl(i+1);
    
    idx_in_interval = find(full_time_history >= t_interval_start & full_time_history <= t_interval_end);
    
    if isempty(idx_in_interval), continue; end
    
    [max_g, idx_max_local] = max(full_g_values(idx_in_interval));
    
    g_max_per_interval(i) = max_g;
    idx_max_global = idx_in_interval(idx_max_local);
    t_max_per_interval(i) = full_time_history(idx_max_global);
end

end