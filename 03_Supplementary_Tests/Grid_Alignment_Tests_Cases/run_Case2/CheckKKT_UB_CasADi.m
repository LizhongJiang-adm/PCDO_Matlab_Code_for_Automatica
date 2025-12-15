function [isConverged, kktValue] = CheckKKT_UB_CasADi(w_opt, grad_f, c_all, grad_c_all, grad_ceq_all, lbw, ubw, strategy, options)
% =========================================================================
% Check KKT Stationarity Conditions for the Upper-Bound CasADi Method
% =========================================================================
% Description:
% Checks the KKT stationarity condition for the full NLP problem transcribed
% by the CasADi-based strategy. It considers active path constraints (at
% their true maximums), other inequality constraints, equality constraints,
% and variable bounds.
%
% Inputs:
%   w_opt        - The optimal decision variable vector.
%   grad_f       - Gradient of the objective function.
%   c_all        - Vector of all inequality constraints.
%   grad_c_all   - Jacobian of all inequality constraints.
%   grad_ceq_all - Jacobian of all equality constraints.
%   lbw, ubw     - Lower and upper bounds for the decision variables.
%   strategy     - The transcription strategy object.
%   options      - The solver options struct (for tolerances).

% Outputs:%
%   isConverged  - A flag (true if converged, false otherwise).
%   kktValue     - The norm of the residual of the KKT stationarity equation.
% =========================================================================

epsilon_act = options.epsilon_act;
stationarityTol = options.epsilon_sat;

%% --- 1. Identify Active Path Constraints at their True Maximums ---
% Get the maximum violation of the true path constraints
[PathCnstrMaxVal, PathCnstrMaxTime_nmlz, ~, ~] = strategy.getMaxPConIntvls(w_opt);

% Find the active time points for each path constraint
ActiveTimeNmlz = cell(size(PathCnstrMaxVal));
for i = 1:numel(PathCnstrMaxVal)
    ActiveTimeNmlz{i} = PathCnstrMaxTime_nmlz{i}(PathCnstrMaxVal{i} >= -epsilon_act);
end

% Get the gradients of the path constraints at these specific active time points
[funPathCnstrDis, funGradPathCnstrDis] = strategy.getPathCnstrDis(ActiveTimeNmlz);
GradPathDis_active = [];
for i = 1:numel(funGradPathCnstrDis)
    if ~isempty(funGradPathCnstrDis{i}) % Only if there are active points for this constraint
        GradPathDis_active = [GradPathDis_active, full(funGradPathCnstrDis{i}(w_opt))'];
    end
end

%% --- 2. Identify Other Active Inequality and Bound Constraints ---
% Find active standard inequality constraints (e.g., terminal constraints)
num_ub_constr = strategy.get_num_ub_constraints();
num_other_ineq = length(c_all) - num_ub_constr;
idx_active_other_ineq = find(c_all(1:num_other_ineq) >= -epsilon_act);
grad_other_ineq_active = grad_c_all(:, idx_active_other_ineq);

% Find active variable bound constraints
idx_active_lb = find(lbw - w_opt >= -epsilon_act);
idx_active_ub = find(w_opt - ubw >= -epsilon_act);

grad_bounds_active = [];
for i = 1:length(idx_active_lb)
    grad = zeros(length(w_opt), 1);
    grad(idx_active_lb(i)) = -1; % For constraint lb - w <= 0
    grad_bounds_active = [grad_bounds_active, grad];
end
for i = 1:length(idx_active_ub)
    grad = zeros(length(w_opt), 1);
    grad(idx_active_ub(i)) = 1; % For constraint w - ub <= 0
    grad_bounds_active = [grad_bounds_active, grad];
end

%% --- 3. Assemble Linear System and Solve for Lagrange Multipliers ---
% Assemble the full Jacobian of active constraints
GradIneq_Active = [grad_other_ineq_active, GradPathDis_active, grad_bounds_active];
GradEq = grad_ceq_all;
C_matrix = [GradIneq_Active, GradEq];
d_vector = -grad_f;

% Set lower bounds for the multipliers (lambda_ineq >= 0)
num_ineq = size(GradIneq_Active, 2);
num_eq = size(GradEq, 2);
lb_lambda = [zeros(num_ineq, 1); -inf(num_eq, 1)];

% Solve the constrained least-squares problem
lsqlin_options = optimoptions('lsqlin', 'Display', 'off');
[lambda, ~, residual_norm_sq] = lsqlin(C_matrix, d_vector, [], [], [], [], lb_lambda, [], [], lsqlin_options);

kktValue = norm(residual_norm_sq); % lsqlin returns squared norm
fprintf('KKT stationarity check: Norm = %.4e (Tolerance = %.4e)\n', kktValue, stationarityTol);

%% --- 4. Determine Convergence ---
if kktValue <= stationarityTol
    isConverged = true;
else
    isConverged = false;
end

end