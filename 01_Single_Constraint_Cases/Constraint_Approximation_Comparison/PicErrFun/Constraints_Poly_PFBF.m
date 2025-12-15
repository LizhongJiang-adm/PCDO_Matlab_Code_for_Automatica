function [c, ceq, grad_c, grad_ceq] = Constraints_Poly_PFBF(decisionVars, problem, T_Dis, ResPmtr)
% =========================================================================
% Polynomial Approximation Constraint Function for PFBF (Exact Replication Version)
% =========================================================================
% This version is a direct, line-by-line replication of the logic in the
% original 'CnstrDisPoly_PFBFc1V5.m' to ensure identical numerical output.

%% --- 1. Parameter Extraction & Initialization (Replicated) ---
N = problem.control.N;
x0 = problem.state.x0;
controlGrid = problem.control.grid;
numStates = problem.system.dimensions;
odeFun = problem.functions.odeSensitivity;

T_Dis = unique(T_Dis(:)');

%% --- 2. Integrate System and Sensitivities ---
% This section is identical to the VDPO version
numSensitivityStates = numStates * N;
y0 = [x0, zeros(1, numSensitivityStates)];
tout = []; xout = [];
opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-7); % Using tolerance from PFBF original

for i = 1:N
    t_start = controlGrid(i); t_end = controlGrid(i+1);
    u_i = decisionVars(i);
    ode_handle = @(t, y) odeFun(t, y, u_i, i, N, [t_start, t_end]);
    [t_interval, x_interval] = ode45(ode_handle, [t_start, t_end], y0, opts);
    y0 = x_interval(end, :);
    if isempty(tout), tout = t_interval; xout = x_interval;
    else, tout = [tout; t_interval(2:end)]; xout = [xout; x_interval(2:end, :)]; end
end

%% --- 3. Compute Pointwise Constraints at T_Dis ---
state_at_T_Dis = interp1(tout, xout, T_Dis, 'PCHIP');
x2_at_T_Dis = state_at_T_Dis(:, 2);

% Pointwise constraints: g(t_i) = x2(t_i) - 0.5 <= 0
c_points = x2_at_T_Dis' - 0.5; % Ensure row vector

% Gradients of pointwise constraints: d(x2)/du
sens_block = state_at_T_Dis(:, numStates+1:end);
grad_c_points = zeros(N, length(T_Dis));
for i = 1:length(T_Dis)
    sens_row = sens_block(i,:);
    x2_q = sens_row(2:numStates:end);
    grad_c_points(:, i) = x2_q';
end

%% --- 4. Compute Polynomial Upper-Bound Constraints (Corrected Logic) ---
poly_intervals = unique([problem.time.T0, T_Dis, problem.time.TF]);
state_at_poly_nodes = interp1(tout, xout, poly_intervals, 'PCHIP')';

num_poly_constraints = length(poly_intervals) - 1;
c_poly = zeros(1, num_poly_constraints);
grad_c_poly = zeros(N, num_poly_constraints);

for i = 1:num_poly_constraints
    t_ends = poly_intervals(i:i+1);
    
    state_vec_left = state_at_poly_nodes(:, i);
    state_vec_right = state_at_poly_nodes(:, i+1);
    
    % g(t) = x2(t) - 0.5
    g_ends = [state_vec_left(2) - 0.5, state_vec_right(2) - 0.5];
    
    % grad(g) = grad(x2)
    sens_row_g_left = state_vec_left(numStates+1:end);
    sens_row_g_right = state_vec_right(numStates+1:end);
    grad_g_ends = [sens_row_g_left(2:numStates:end), sens_row_g_right(2:numStates:end)];
    
    % --- Calculate the correct one-sided derivatives and their gradients ---
    dot_g_ends = zeros(1, 2);
    grad_dot_g_ends = zeros(N, 2);
    
    % --- Info for LEFT endpoint t_ends(1): Use RIGHT-SIDED derivative g_dot(t+) ---
    control_idx_left = find(controlGrid <= t_ends(1) + 1e-10, 1, 'last');
    if t_ends(1) == controlGrid(end), control_idx_left = N; end
    u_node_left = decisionVars(control_idx_left);
    dxdt_vec_left = odeFun(t_ends(1), state_vec_left, u_node_left, control_idx_left, N, []);
    dot_g_ends(1) = dxdt_vec_left(2); % g_dot = x2_dot
    sens_dxdt_left = dxdt_vec_left(numStates+1:end);
    grad_dot_g_ends(:, 1) = sens_dxdt_left(2:numStates:end);

    % --- Info for RIGHT endpoint t_ends(2): Use LEFT-SIDED derivative g_dot(t-) ---
    control_idx_right = find(controlGrid < t_ends(2) - 1e-10, 1, 'last');
    if isempty(control_idx_right), control_idx_right = 1; end
    u_node_right = decisionVars(control_idx_right);
    dxdt_vec_right = odeFun(t_ends(2), state_vec_right, u_node_right, control_idx_right, N, []);
    dot_g_ends(2) = dxdt_vec_right(2); % g_dot = x2_dot
    sens_dxdt_right = dxdt_vec_right(numStates+1:end);
    grad_dot_g_ends(:, 2) = sens_dxdt_right(2:numStates:end);
    
    [c_poly(i), grad_c_poly(:, i)] = PloyMax_Gra(t_ends, g_ends, dot_g_ends, grad_g_ends, grad_dot_g_ends);
end

c_poly = c_poly + reshape(ResPmtr,1,[]);

%% --- 5. Combine All Constraints ---
c = [c_points, c_poly];
grad_c = [grad_c_points, grad_c_poly];
ceq = [];
grad_ceq = [];

end