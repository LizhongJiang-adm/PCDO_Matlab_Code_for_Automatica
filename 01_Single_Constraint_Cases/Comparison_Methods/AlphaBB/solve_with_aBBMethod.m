function [solution, stats] = solve_with_aBBMethod(problem, options)
% =========================================================================
% Solve OCP using the Iterative aBB Method
% =========================================================================
% Description:
% This function implements an iterative algorithm for solving optimal
% control problems using the aBB method. In each iteration, it dynamically
% constructs an extended NLP subproblem that includes control variables (u),
% time variables (t_k), and slack variables (gamma). It refines a time
% grid ('eta') until KKT conditions are met.
%
% Inputs:
%   problem - A struct containing the problem definition (control part).
%   options - A struct containing solver settings and tolerances.
%
% Outputs:
%   solution - A struct containing the final optimization results.
%   stats    - A struct with statistics about the solution process.
% =========================================================================

%% --- 1. Initialization ---
startTime = tic;
iterationCount = 1;

% Extract key parameters from inputs
N = problem.numVars;
% objFun_handle = @(u, cGrid) problem.functions.objective(u(1:N), problem.control.grid, problem.state.x0);
constraintsFun = @(u, cGrid) problem.functions.constraints_aBB(u, problem, cGrid);
disConstrFun = @(u, tMax) problem.functions.constraintsDis(u, problem.control.grid, problem.state.x0, tMax);
objFun_handle = @(u_comp) objective_wrapper_aBB(u_comp, problem);
% To ensure the final solution strictly satisfies g(u) <= 0 rather than g(u) <= Tol,
% we shift the constraints passed to fmincon by its own tolerance.
fminconConstrFun = @(u, cGrid) shifted_constraint_wrapper(u, cGrid, problem, options.fmincon.ConstraintTolerance);


nlpOptions = options.fmincon;

% --- Algorithm State Initialization ---
% The state of the aBB algorithm includes the control variables and the grid.
current_control_vars = problem.initialGuess_control;
eta = problem.control.grid; % Start with the control grid as the initial eta

% History tracking
history.eta{iterationCount} = eta;
history.controlVars{iterationCount} = current_control_vars;

fprintf('\n--- Starting aBB Method Solver ---\n');

%% --- 2. Main Iteration Loop ---
while true
    fprintf('\n--- Iteration: %d ---\n', iterationCount);
    fprintf('Number of eta grid points: %d\n', length(eta));

    % ================================================================
    % --- Dynamically Construct the NLP Subproblem for the current eta ---
    % This logic comes from the original 'subproblem_aBB.m'.
    % ================================================================
    M_eta = length(eta) - 1;
    num_total_vars = N + 3*M_eta;
    
    % 1. Build the composite initial guess (u0)
    u0 = zeros(num_total_vars, 1);
    u0(1:N) = current_control_vars;
    u0(N+1 : N+M_eta) = 0.5 * (eta(1:end-1)' + eta(2:end)'); % Midpoints for t_k
    u0(N+M_eta+1 : end) = 0.5; % Initial guess for gammas
    
    % 2. Build the composite lower and upper bounds (LB, UB)
    LB = zeros(num_total_vars, 1);
    UB = zeros(num_total_vars, 1);
    LB(1:N) = problem.lowerBounds_control;
    UB(1:N) = problem.upperBounds_control;
    LB(N+1 : N+M_eta) = eta(1:end-1)';
    UB(N+1 : N+M_eta) = eta(2:end)';
    LB(N+M_eta+1 : end) = -Inf;
    UB(N+M_eta+1 : end) = Inf;

    % 3. Create function handles for fmincon, passing the current eta
    % obj_handle = @(u_comp) objFun(u_comp(1:N), problem);
    % constr_handle = @(u_comp) constraintsFun(u_comp, problem, eta);

    % --- Solve the NLP Subproblem ---
    [decvar, objVal, exitFlag, output, lambda, grad] = fmincon(...
        @(u)objFun_handle(u), u0, [], [], [], [], LB, UB, @(u)fminconConstrFun(u,eta), nlpOptions);
    
    history.cpuTime(iterationCount) = toc(startTime);
    
    if exitFlag > 0
        fprintf('NLP subproblem solved successfully.\n');
        % Always update the starting point for the next iteration (warm start)
        current_control_vars = decvar(1:N);
        history.isFeasible(iterationCount) = 1;
    else
        fprintf('Warning: NLP subproblem failed or was infeasible (Exit Flag: %d).\n', exitFlag);
        history.isFeasible(iterationCount) = 0;
    end
    
    % --- Convergence Check and Grid Refinement ---
    
    % Step A: Check KKT conditions
    grad_f_control = grad(1:N); % Gradient of objective w.r.t control vars
    decvar_control = decvar(1:N); 
    % --- Check for Convergence ---
    t_maxpoints_cell = {decvar(N+1 : N+M_eta)};
    % Step B: Evaluate the discrete path constraints at these midpoints.
    [cieq_Dis, ~, gCieq_Dis, ~] = disConstrFun(decvar_control, t_maxpoints_cell);
    % Step C: Identify active constraints based on their values at the midpoints.
    activePathIdx = find(cieq_Dis > -options.activeConstraintTol & cieq_Dis <= 0);
    gCieqActive = gCieq_Dis(:, activePathIdx);
    [isConverged, kktValue] = KKT_MinStat_aBB(decvar_control, grad_f_control, gCieqActive, problem, options);

    
    % [isConverged, kktInfo] = CheckKKT_aBB(decvar, grad_f_control, problem, options); % Placeholder
    
    if isConverged && exitFlag > 0 % Converged and feasible
        fprintf('\nConvergence criteria met. Terminating.\n');
        break;
    end
    
    % Step B: If not converged, identify intervals to refine
    if exitFlag < 0
        % If NLP is infeasible, refine all intervals
        indicesToRefine = 1:M_eta;
    end
    
    % Step 1: Identify which intervals to refine based on active constraints.
    % This logic comes from your original main script.
    [c, ~, ~, ~] = fminconConstrFun(decvar,eta); % c are the凹化constraint values
    indicesToRefine = find(c > - options.activeConstraintTol)'; % Ensure row vector
    
    % Step 2: Extract the corresponding time points 't_k' for these active intervals.
    active_tk = decvar(N + indicesToRefine);
    
    % Step 3: Call the final, definitive refinement function.
    % Note that it now takes 'options' and 'problem' as input.
    [new_eta, control_flow] = refine_grid_aBB(eta, indicesToRefine, active_tk, exitFlag, options, problem);
    
    if control_flow == 1
        fprintf('Refinement function signaled to terminate the algorithm due to stagnation.\n');
        break;
    elseif control_flow == 2
        fprintf('Refinement function signaled an unexpected error. Terminating.\n');
        break;
    end

    % --- Update State for Next Iteration ---
    eta = new_eta;
    iterationCount = iterationCount + 1;
    history.eta{iterationCount} = eta;
    history.controlVars{iterationCount} = current_control_vars;
    
    if iterationCount > 100 % Safety break
        fprintf('Warning: Maximum number of iterations (100) reached. Terminating.\n');
        break;
    end
end

%% --- 3. Package Results ---
fprintf('\n--- Algorithm Finished ---\n');
totalTime = toc(startTime);

solution.decisionVariables = decvar(1:N);
solution.timeVariables_tk = decvar(N+1 : N+M_eta);
solution.gamma_L = decvar(N+M_eta+1 : N+2*M_eta);
solution.gamma_U = decvar(N+2*M_eta+1 : N+3*M_eta);
solution.fullDecisionVector = decvar;
solution.objectiveValue = objVal;
solution.exitFlag = exitFlag;

stats.totalIterations = iterationCount;
stats.cpuTime = totalTime;
stats.fminconOutput = output;
stats.history = history;

fprintf('Total iterations: %d\n', stats.totalIterations);
fprintf('Total CPU time: %.2f seconds\n', stats.cpuTime);

end

function [J, grad_J_full] = objective_wrapper_aBB(u_composite, problem)
% =========================================================================
% Wrapper for the aBB objective function to handle extended gradients.
% =========================================================================
% Description:
% This function calls the original objective function, which only depends on
% the first N control variables. It then pads the resulting gradient vector
% with zeros to match the full dimension of the aBB method's decision
% vector [u_control; t_k; gamma_L; gamma_U].
%
% Inputs:
%   u_composite - The full decision vector (N + 3*M_eta).
%   problem     - The problem definition struct.
%
% Outputs:
%   J           - The scalar objective value.
%   grad_J_full - The full gradient vector, padded with zeros.
% =========================================================================

    % Extract the number of control variables
    N = problem.numVars;
    
    % Extract only the control variables to pass to the original objective function
    u_control = u_composite(1:N);
    
    % Call the original objective function to get the value and the partial gradient
    [J, grad_J_control] = problem.functions.objective(u_control, problem.control.grid, problem.state.x0);
    
    % Get the total number of variables in the composite vector
    num_total_vars = length(u_composite);
    
    % Initialize the full gradient vector with zeros
    grad_J_full = zeros(num_total_vars, 1);
    
    % Place the non-zero part of the gradient in the first N positions
    % Ensure the gradient has the same column/row orientation as the zero vector
    grad_J_full(1:N) = grad_J_control(:);

end

function [c, ceq, grad_c, grad_ceq] = shifted_constraint_wrapper(u, cGrid, problem, shift_value)
    % Call the original upper-bound constraint function
    [c_original, ceq, grad_c, grad_ceq] = problem.functions.constraints_aBB(u, problem, cGrid);
    
    % Apply the upward shift to the inequality constraint values.
    % The gradient remains unchanged because the derivative of a constant (shift_value) is zero.
    c = c_original + shift_value;
end