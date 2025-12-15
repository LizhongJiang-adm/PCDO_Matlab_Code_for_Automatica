function t_max_poly = GetPolyTmax_PerInterval_Poly(decisionVars, problem, T_Dis)
% =========================================================================
% Find Time of Maximum for Polynomial Approximations (Generic, Corrected Version)
% =========================================================================
% Description:
% For each interval defined by nodes [t0, T_Dis, tf], this function
% constructs a Hermite cubic polynomial approximation of the path constraint.
% It then finds and returns the time 't_max' at which this polynomial
% reaches its maximum within each interval. This version correctly handles
% derivative discontinuities at control switching points and is fully generic.
%
% Inputs:
%   decisionVars - The vector of control variables (u).
%   problem      - The problem definition struct.
%   T_Dis        - A vector of discrete time points.
%
% Outputs:
%   t_max_poly   - A vector containing the time of the polynomial's
%                  maximum for each interval.
% =========================================================================

%% --- 1. Extract Parameters & Initialize ---
N = problem.control.N;
x0 = problem.state.x0;
controlGrid = problem.control.grid;
numStates = problem.system.dimensions;
odeFun = problem.functions.odeSensitivity;
pathConstraintHandle = problem.pathConstraints.handle{1};
dotPathConstraintHandle = problem.functions.dotPathConstraint;

%% --- 2. Integrate System and Sensitivities ---
% Perform one single, high-fidelity integration run to get the full trajectory.
numSensitivityStates = numStates * N;
y0 = [x0, zeros(1, numSensitivityStates)];
tout = []; xout = [];
opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

for i = 1:N
    t_start = controlGrid(i); t_end = controlGrid(i+1);
    u_i = decisionVars(i);
    ode_handle = @(t, y) odeFun(t, y, u_i, i, N, [t_start, t_end]);
    [t_interval, x_interval] = ode45(ode_handle, [t_start, t_end], y0, opts);
    y0 = x_interval(end, :);
    if isempty(tout), tout = t_interval; xout = x_interval;
    else, tout = [tout; t_interval(2:end)]; xout = [xout; x_interval(2:end, :)]; end
end

%% --- 3. Construct Polynomials and Find t_max (Corrected & Generic Logic) ---
poly_intervals = unique([problem.time.T0, T_Dis(:)', problem.time.TF]);
state_at_poly_nodes = interp1(tout, xout, poly_intervals, 'PCHIP')';

num_poly_intervals = length(poly_intervals) - 1;
t_max_poly = zeros(1, num_poly_intervals);

for i = 1:num_poly_intervals
    t_ends = poly_intervals(i:i+1);
    
    % Get state vectors at the interval endpoints
    state_vec_left = state_at_poly_nodes(:, i);
    state_vec_right = state_at_poly_nodes(:, i+1);
    
    dot_g_ends = zeros(1, 2);
    
    % --- Calculate info for LEFT endpoint t_ends(1): Use RIGHT-SIDED derivative g_dot(t+) ---
    control_idx_left = find(controlGrid <= t_ends(1) + 1e-10, 1, 'last');
    if t_ends(1) == controlGrid(end), control_idx_left = N; end
    u_node_left = decisionVars(control_idx_left);
    
    % --- Calculate info for RIGHT endpoint t_ends(2): Use LEFT-SIDED derivative g_dot(t-) ---
    control_idx_right = find(controlGrid < t_ends(2) - 1e-10, 1, 'last');
    if isempty(control_idx_right), control_idx_right = 1; end
    u_node_right = decisionVars(control_idx_right);

    % --- GENERIC CALCULATION of g and g_dot using function handles ---
    g_ends = [pathConstraintHandle(t_ends(1), state_vec_left(1:numStates)', u_node_left), ...
              pathConstraintHandle(t_ends(2), state_vec_right(1:numStates)', u_node_right)];
              
    dxdt_vec_left = odeFun(t_ends(1), state_vec_left, u_node_left, control_idx_left, N, []);
    dot_g_ends(1) = dotPathConstraintHandle(t_ends(1), dxdt_vec_left');
    
    dxdt_vec_right = odeFun(t_ends(2), state_vec_right, u_node_right, control_idx_right, N, []);
    dot_g_ends(2) = dotPathConstraintHandle(t_ends(2), dxdt_vec_right');
    
    % Call PloyMax helper to get tmax
    [~,~,~,~,tmax,~] = PloyMax(t_ends, g_ends, dot_g_ends);

    % Clamp tmax to be within the interval boundaries
    if tmax <= t_ends(1)
        t_max_poly(i) = t_ends(1);
    elseif tmax >= t_ends(2)
        t_max_poly(i) = t_ends(2);
    else
        t_max_poly(i) = tmax;
    end
end

end