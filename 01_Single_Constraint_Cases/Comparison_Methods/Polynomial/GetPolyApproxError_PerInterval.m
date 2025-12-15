function [poly_error_max, t_at_max_error] = GetPolyApproxError_PerInterval(decisionVars, problem, T_Dis)
% =========================================================================
% Find Max Error Between Path Constraint and Polynomial Approximation
% =========================================================================
% Description:
% For each interval defined by the nodes [t0, T_Dis, tf], this function
% computes the true path constraint g(t) and its Hermite cubic polynomial
% approximation p(t). It then finds the maximum of the absolute error
% |g(t) - p(t)| and the time at which it occurs for each interval.
%
% Inputs:
%   decisionVars - The vector of control variables (u).
%   problem      - The problem definition struct.
%   T_Dis        - A vector of discrete time points defining the intervals.
%
% Outputs:
%   poly_error_max   - A vector of the max error for each interval.
%   t_at_max_error   - A vector of the time of max error for each interval.
% =========================================================================

%% --- 1. Extract Parameters & Initialize ---
N = problem.control.N;
x0 = problem.state.x0;
controlGrid = problem.control.grid;
numStates = problem.system.dimensions;
odeFun = problem.functions.odeSensitivity;
% For g_dot calculation, we need a state-only ODE function.
% Let's assume the problem definition provides this.
odeStateFun = @(t, y, u, idx) odeFun(t, [y; zeros(numStates*N, 1)], u, idx, N, []);

pathConstraintHandle = problem.pathConstraints.handle{1};
dotPathConstraintHandle = problem.functions.dotPathConstraint;

%% --- 2. Simulate Trajectory (with Sensitivities) ---
% We perform one single, high-fidelity integration run.
numSensitivityStates = numStates * N;
y0 = [x0, zeros(1, numSensitivityStates)];
tout = []; xout = [];
% opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);, opts

for i = 1:N
    t_start = controlGrid(i); t_end = controlGrid(i+1);
    u_i = decisionVars(i);
    ode_handle = @(t, y) odeFun(t, y, u_i, i, N, [t_start, t_end]);
    [t_interval, x_interval] = ode45(ode_handle, [t_start, t_end], y0);
    y0 = x_interval(end, :);
    if isempty(tout)
        tout = t_interval; xout = x_interval;
    else
        tout = [tout; t_interval(2:end)]; xout = [xout; x_interval(2:end, :)];
    end
end

%% --- 3. Analyze Each Interval ---
poly_intervals = unique([problem.time.T0, T_Dis(:)', problem.time.TF]);
num_poly_intervals = length(poly_intervals) - 1;
poly_error_max = zeros(1, num_poly_intervals);
t_at_max_error = zeros(1, num_poly_intervals);

for i = 1:num_poly_intervals
    t_interval_start = poly_intervals(i);
    t_interval_end = poly_intervals(i+1);
    
    % --- A. Get True Trajectory in the Interval ---
    idx_in_interval = find(tout >= t_interval_start & tout <= t_interval_end);
    if isempty(idx_in_interval), continue; end
    
    t_subset = tout(idx_in_interval);
    state_subset = xout(idx_in_interval, 1:numStates);
    
    % Reconstruct control for this subset
    u_subset = zeros(size(t_subset));
    for k = 1:length(t_subset)
        control_idx_k = find(controlGrid <= t_subset(k) + 1e-10, 1, 'last');
        if t_subset(k) == controlGrid(end), control_idx_k = N; end
        u_subset(k) = decisionVars(control_idx_k);
    end
    
    % Evaluate true path constraint g(t) in the interval
    g_true_values = pathConstraintHandle(t_subset, state_subset, u_subset);

    % --- B. Construct Polynomial Approximation p(t) (CORRECTED) ---
    t_ends = [t_interval_start, t_interval_end];
    
    % Get state at endpoints
    state_ends_matrix = interp1(tout, xout, t_ends, 'PCHIP');
    state_ends = [state_ends_matrix(1, 1:numStates)', state_ends_matrix(2, 1:numStates)'];
    
    g_ends = [pathConstraintHandle(t_ends(1), state_ends(:,1)', NaN), ...
              pathConstraintHandle(t_ends(2), state_ends(:,2)', NaN)];

    % =========================== KEY CORRECTION HERE ===========================
    % Calculate the correct one-sided derivatives for each endpoint.
    dot_g_ends = zeros(1, 2);
    
    % --- Derivative at LEFT endpoint t_ends(1): Use RIGHT-SIDED derivative g_dot(t+) ---
    control_idx_left = find(controlGrid <= t_ends(1) + 1e-10, 1, 'last');
    if t_ends(1) == controlGrid(end), control_idx_left = N; end
    u_node_left = decisionVars(control_idx_left);
    dxdt_vec_left = odeStateFun(t_ends(1), state_ends(:,1), u_node_left, control_idx_left);
    dot_g_ends(1) = dotPathConstraintHandle(t_ends(1), dxdt_vec_left');

    % --- Derivative at RIGHT endpoint t_ends(2): Use LEFT-SIDED derivative g_dot(t-) ---
    control_idx_right = find(controlGrid < t_ends(2) - 1e-10, 1, 'last');
    if isempty(control_idx_right), control_idx_right = 1; end
    u_node_right = decisionVars(control_idx_right);
    dxdt_vec_right = odeStateFun(t_ends(2), state_ends(:,2), u_node_right, control_idx_right);
    dot_g_ends(2) = dotPathConstraintHandle(t_ends(2), dxdt_vec_right');
    % ===========================================================================
    
    g_poly_values = GetPloyVal(t_ends, g_ends, dot_g_ends, t_subset);

    % --- C. Find Max Error ---
    abs_error = abs(g_true_values - g_poly_values);
    [max_err, idx_max_local] = max(abs_error);
    
    poly_error_max(i) = max_err;
    t_at_max_error(i) = t_subset(idx_max_local);
end

end