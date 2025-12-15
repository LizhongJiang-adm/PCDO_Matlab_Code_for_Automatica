function [c, ceq, grad_c, grad_ceq] = Constraints_Intvl_JCB(decisionVars, constraintGrid, problem)
% =========================================================================
% Interval Method Constraint Function for JCB Problem (Standardized)
% =========================================================================
% Description:
% Computes the nonlinear inequality constraints (c) and their gradients (grad_c)
% for the NLP subproblem in the Interval Method for the JCB problem.
% This function follows the standardized, self-contained format.
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
numStates = problem.system.dimensions; % sys_dms for JCB
UBdot2g = problem.algorithm_params.UBdot2g;
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
    if i == numVars
        active_constraint_indices = find(intervalMidpoints >= t_control_start & intervalMidpoints <= t_control_end);
    end

    if isempty(active_constraint_indices)
        [~, x_out] = ode45(@(t,x) odeFun(t, x, controlParams(i), i, numVars, [t_control_start, t_control_end]), ...
                         [t_control_start, t_control_end], y0, opts);
        y0 = x_out(end,:);
        continue;
    end

    [t_ode, x_ode] = ode45(@(t,x) odeFun(t, x, controlParams(i), i, numVars, [t_control_start, t_control_end]), ...
                         [t_control_start, t_control_end], y0, opts);
    y0 = x_ode(end,:);

    midpoints_in_interval = intervalMidpoints(active_constraint_indices);
    state_at_midpoints = interp1(t_ode, x_ode, midpoints_in_interval, 'PCHIP')';

    % --- Evaluate Constraints and Gradients at Midpoints for JCB PROBLEM ---
    % Path constraint: g(t,x) = x2 + 0.5 - 8*(t - 0.5)^2
    % Time derivative: g_dot(t,x,u) = u - x2 - 16*(t - 0.5)
    
    x2 = state_at_midpoints(2,:);
    
    g0 = x2 + 0.5 - 8 .* (midpoints_in_interval - 0.5).^2;
    Nabla_g0 = controlParams(i) - x2 - 16 .* (midpoints_in_interval - 0.5);
    
    w_subT_active = intervalWidths(active_constraint_indices);
    c(active_constraint_indices) = g0 + 0.5 * Smoothing_absolute(Nabla_g0) .* w_subT_active ...
        + 0.5 * UBdot2g * 0.25 * w_subT_active.^2;
        
    % --- Calculate Gradient of Constraints ---
    g0_u_grad = zeros(numVars, length(active_constraint_indices));
    Nabla_g0_u_grad = zeros(numVars, length(active_constraint_indices));

    for p = 1:numVars % Gradient w.r.t. each u_p
        % Extract sensitivity of x2 w.r.t. u_p
        % For JCB, x is [x1, x2]. So numStates=2. Sensitivity of x2 is 2, 2+2, 2+2+2, ...
        x2_q_p = state_at_midpoints(numStates*p + 2, :);
        
        % Gradient of g0 = x2_q_p
        g0_u_grad(p, :) = x2_q_p;
        
        % Gradient of Nabla_g0 = -x2_q_p
        Nabla_g0_u_grad(p, :) = -x2_q_p;
        
        if p == i
            % Add the direct derivative w.r.t. u_i
            Nabla_g0_u_grad(p, :) = Nabla_g0_u_grad(p, :) + 1;
        end
    end
    
    Sms_Nabla_g0_u = 0.5 * w_subT_active .* Smoothing_absolute_der(Nabla_g0) .* Nabla_g0_u_grad;
    
    grad_c(:, active_constraint_indices) = g0_u_grad + Sms_Nabla_g0_u;
end

end