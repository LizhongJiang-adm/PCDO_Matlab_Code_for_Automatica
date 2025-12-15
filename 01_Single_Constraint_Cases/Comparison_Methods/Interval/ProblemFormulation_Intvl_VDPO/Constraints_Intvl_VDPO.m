function [c, ceq, grad_c, grad_ceq] = Constraints_Intvl_VDPO(decisionVars, constraintGrid, problem)
% =========================================================================
% Interval Method Constraint Function for VDPO (Refactored)
% =========================================================================
% Description:
% Computes the nonlinear inequality constraints (c) and their gradients (grad_c)
% for the NLP subproblem in the Interval Method. The constraints are built
% upon the grid provided by 'constraintGrid'. This version is self-contained.
%
% Inputs:
%   decisionVars   - The vector of decision variables (u).
%   constraintGrid - The current time grid for constructing constraints (subdvd_T).
%   problem        - The problem definition struct.
%
% Outputs:
%   c, ceq         - Inequality and equality constraint vectors.
%   grad_c, grad_ceq - Gradients of the constraints.
% =========================================================================

%% --- 1. Extract Parameters from Inputs ---
x0 = problem.state.x0;
numVars = problem.numVars; % N
controlGrid = problem.control.grid; % CvpGrd
numStates = problem.system.dimensions; % sys_dms
UBdot2g = problem.algorithm_params.UBdot2g;

% Assume the sensitivity ODE function is passed via the problem struct
% It seems a different ODE function is used here compared to the objective
% Let's add it to the problem definition later.
% For now, let's assume it's problem.functions.odeSensitivity.
odeFun = problem.functions.odeSensitivity;

numConstraintIntervals = length(constraintGrid) - 1;

%% --- 2. Initialize Outputs ---
c = zeros(1, numConstraintIntervals);
ceq = [];
grad_c = zeros(numVars, numConstraintIntervals);
grad_ceq = [];

%% --- 3. Main Calculation Loop ---
% Initialize the full state vector for integration [states, sensitivities]
numSensitivityStates = numStates * numVars;
y0 = [x0, zeros(1, numSensitivityStates)];

% Pre-calculate interval midpoints and widths for efficiency
intervalMidpoints = 0.5 * (constraintGrid(1:end-1) + constraintGrid(2:end));
intervalWidths = constraintGrid(2:end) - constraintGrid(1:end-1);

opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

% Reshape decision vars for easier access
controlParams = reshape(decisionVars, 1, numVars);

% Loop over each CONTROL interval to perform integration
for i = 1:numVars % Loop from 1 to N
    t_control_start = controlGrid(i);
    t_control_end = controlGrid(i+1);
    
    % Find which constraint interval midpoints fall within this control interval
    active_constraint_indices = find(intervalMidpoints >= t_control_start & intervalMidpoints < t_control_end);
    % Handle the endpoint case
    if i == numVars
        active_constraint_indices = find(intervalMidpoints >= t_control_start & intervalMidpoints <= t_control_end);
    end

    if isempty(active_constraint_indices)
        % If no midpoints are in this control interval, we still need to integrate
        % to get the correct initial state for the next interval.
        [~, x_out] = ode45(@(t,x) odeFun(t, x, controlParams(i), i, numVars, [t_control_start, t_control_end]), ...
                         [t_control_start, t_control_end], y0, opts);
        y0 = x_out(end,:);
        continue; % Move to the next control interval
    end

    % Integrate across the current control interval
    [t_ode, x_ode] = ode45(@(t,x) odeFun(t, x, controlParams(i), i, numVars, [t_control_start, t_control_end]), ...
                         [t_control_start, t_control_end], y0, opts);
    y0 = x_ode(end,:);

    % Interpolate to get state and sensitivity values at the midpoints
    midpoints_in_interval = intervalMidpoints(active_constraint_indices);
    state_at_midpoints = interp1(t_ode, x_ode, midpoints_in_interval, 'PCHIP')';

    % --- Evaluate Constraints and Gradients at Midpoints ---
    % Path constraint: g(x) = -x1 - 0.4
    % Time derivative of path constraint: g_dot(x,u) = -u + x2 + x1*(x2^2 - 1)
    
    x1 = state_at_midpoints(1,:);
    x2 = state_at_midpoints(2,:);
    
    g0 = -x1 - 0.4; % Value of g at midpoints
    Nabla_g0 = -controlParams(i) + x2 + x1 .* (x2.^2 - 1); % Value of g_dot at midpoints
    
    % Calculate the constraint value using the interval formula
    w_subT_active = intervalWidths(active_constraint_indices);
    c(active_constraint_indices) = g0 + 0.5 * Smoothing_absolute(Nabla_g0) .* w_subT_active ...
        + 0.5 * UBdot2g * 0.25 * w_subT_active.^2;
        
    % --- Calculate Gradient of Constraints ---
    g0_u_grad = zeros(numVars, length(active_constraint_indices));
    Nabla_g0_u_grad = zeros(numVars, length(active_constraint_indices));

    for p = 1:numVars % Gradient w.r.t. each u_p
        % Extract sensitivities of x1 and x2 w.r.t. u_p
        x1_q_p = state_at_midpoints(numStates*p + 1, :);
        x2_q_p = state_at_midpoints(numStates*p + 2, :);
        
        % Gradient of g0 = -x1_q_p
        g0_u_grad(p, :) = -x1_q_p;
        
        % Gradient of Nabla_g0
        dg_dx1 = x2.^2 - 1;
        dg_dx2 = 1 + 2*x1.*x2;
        Nabla_g0_u_grad(p, :) = dg_dx1 .* x1_q_p + dg_dx2 .* x2_q_p;
        
        if p == i
            % Add the direct derivative w.r.t. u_i
            Nabla_g0_u_grad(p, :) = Nabla_g0_u_grad(p, :) - 1;
        end
    end
    
    Sms_Nabla_g0_u = 0.5 * w_subT_active .* Smoothing_absolute_der(Nabla_g0) .* Nabla_g0_u_grad;
    
    grad_c(:, active_constraint_indices) = g0_u_grad + Sms_Nabla_g0_u;
end

end