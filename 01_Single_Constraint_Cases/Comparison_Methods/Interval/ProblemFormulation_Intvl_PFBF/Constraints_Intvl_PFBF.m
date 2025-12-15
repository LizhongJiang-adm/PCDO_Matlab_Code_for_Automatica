function [c, ceq, grad_c, grad_ceq] = Constraints_Intvl_PFBF(decisionVars, constraintGrid, problem)
% =========================================================================
% Interval Method Constraint Function for PFBF Problem (Standardized)
% =========================================================================
% Description:
% Computes the nonlinear inequality constraints (c) and their gradients (grad_c)
% for the NLP subproblem in the Interval Method for the PFBF problem.
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
numStates = problem.system.dimensions; % sys_dms for PFBF
UBdot2g = problem.algorithm_params.UBdot2g;
odeFun = problem.functions.odeSensitivity;

% Physical parameters for the PFBF model
kL=0.006; mu=0.11; Yxs=0.47; theta=0.004; Yp=1.2;
kI=0.1; Mx=0.029; Sl=400; kXP=0.01; kP=0.0001;

numConstraintIntervals = length(constraintGrid) - 1;

%% --- 2. Initialize Outputs ---
c = zeros(1, numConstraintIntervals);
ceq = [];
grad_c = zeros(numVars, numConstraintIntervals);
grad_ceq = [];

%% --- 3. Main Calculation Loop ---
numSensitivityStates = numStates * numVars;
y0 = [x0, zeros(1, numSensitivityStates)];

intervalMidpoints = 0.5 * (constraintGrid(1:end-1) + constraintGrid(2:end));
intervalWidths = constraintGrid(2:end) - constraintGrid(1:end-1);

opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
controlParams = reshape(decisionVars, 1, numVars);

for i = 1:numVars % Loop over each CONTROL interval
    t_control_start = controlGrid(i);
    t_control_end = controlGrid(i+1);
    
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

    % --- Evaluate Constraints and Gradients at Midpoints for PFBF PROBLEM ---
    x1 = state_at_midpoints(1,:);
    x2 = state_at_midpoints(2,:);
    x3 = state_at_midpoints(3,:);
    x4 = state_at_midpoints(4,:);
    u_i = controlParams(i);
    
    % Path constraint: g(x) = x2 - 0.5
    g0 = x2 - 0.5;
    
    % Time derivative of path constraint: g_dot(x,u)
    term1 = (u_i .* (Sl - x2)) ./ x4;
    term2 = Mx .* x1;
    term3_num = theta .* x1 .* x2;
    term3_den = Yp .* (x2.^2./kI + x2 + kP);
    term4_num = mu .* x1 .* x2;
    term4_den = Yxs .* (x2 + kL .* x1);
    Nabla_g0 = term1 - term2 - term3_num ./ term3_den - term4_num ./ term4_den;
    
    w_subT_active = intervalWidths(active_constraint_indices);
    c(active_constraint_indices) = g0 + 0.5 * Smoothing_absolute(Nabla_g0) .* w_subT_active ...
        + 0.5 * UBdot2g * 0.25 * w_subT_active.^2;
        
    % --- Calculate Gradient of Constraints ---
    g0_u_grad = zeros(numVars, length(active_constraint_indices));
    Nabla_g0_u_grad = zeros(numVars, length(active_constraint_indices));

    for p = 1:numVars % Gradient w.r.t. each u_p
        x1_q_p = state_at_midpoints(numStates*p + 1, :);
        x2_q_p = state_at_midpoints(numStates*p + 2, :);
        x3_q_p = state_at_midpoints(numStates*p + 3, :);
        x4_q_p = state_at_midpoints(numStates*p + 4, :);
        
        % Gradient of g0 = x2_q_p
        g0_u_grad(p, :) = x2_q_p;
        
        % Gradient of Nabla_g0 (term by term)
        dg_dx1 = (kL.*mu.*x1.*x2)./(Yxs.*(x2 + kL.*x1).^2) - (mu.*x2)./(Yxs.*(x2 + kL.*x1)) - (theta.*x2)./(Yp.*(x2.^2./kI + x2 + kP)) - Mx;
        dg_dx2 = (mu.*x1.*x2)./(Yxs.*(x2 + kL.*x1).^2) - (mu.*x1)./(Yxs.*(x2 + kL.*x1)) - (theta.*x1)./(Yp.*(x2.^2./kI + x2 + kP)) - u_i./x4 + (theta.*x1.*x2.*((2.*x2)./kI + 1))./(Yp.*(x2.^2./kI + x2 + kP).^2);
        dg_dx3 = 0;
        dg_dx4 = -(u_i.*(Sl - x2))./x4.^2;

        Nabla_g0_u_grad(p, :) = dg_dx1.*x1_q_p + dg_dx2.*x2_q_p + dg_dx3.*x3_q_p + dg_dx4.*x4_q_p;
        
        if p == i
            % Add the direct derivative w.r.t. u_i
            dg_du = (Sl - x2)./x4;
            Nabla_g0_u_grad(p, :) = Nabla_g0_u_grad(p, :) + dg_du;
        end
    end
    
    Sms_Nabla_g0_u = 0.5 * w_subT_active .* Smoothing_absolute_der(Nabla_g0) .* Nabla_g0_u_grad;
    
    grad_c(:, active_constraint_indices) = g0_u_grad + Sms_Nabla_g0_u;
end

end