function [solution, stats] = solve_with_IntervalMethod(problem, options)
% =========================================================================
% Solve OCP using the Iterative Interval Method
% =========================================================================
% Description:
% This function implements an iterative algorithm for solving optimal
% control problems with path constraints using an interval-based method.
% It calls a standard NLP solver (fmincon) in each iteration and refines
% a time grid ('subdvd_T') until KKT conditions are met.
%
% Inputs:
%   problem - A struct containing the full problem definition.
%   options - A struct containing solver settings and tolerances.
%
% Outputs:
%   solution - A struct containing the final optimization results.
%   stats    - A struct with statistics about the solution process.
% =========================================================================

%% --- 1. Initialization ---

% Start timer
startTime = tic;
iterationCount = 1;

% Extract parameters from input structs to make code cleaner
numVars = problem.numVars;
lowerBounds = problem.lowerBounds;
upperBounds = problem.upperBounds;
initialGuess = problem.initialGuess;
controlGrid = problem.control.grid; % CvpGrd

% Initialize the dynamic state of the solver
current_decision_vars = initialGuess;
% This is the grid that gets refined. Corresponds to 'subdvd_T'.
constraintGrid = controlGrid; 

% Extract function handles
objFun_handle = @(u) problem.functions.objective(u, problem.control.grid, problem.state.x0);
constrFun_handle = @(u, cGrid) problem.functions.constraintsIntvl(u, cGrid, problem);
disConstrFun = @(u, tMax) problem.functions.constraintsDis(u, problem.control.grid, problem.state.x0, tMax);
% To ensure the final solution strictly satisfies g(u) <= 0 rather than g(u) <= Tol,
% we shift the constraints passed to fmincon by its own tolerance.
fminconConstrFun = @(u, cGrid) shifted_constraint_wrapper(u, cGrid, problem, options.fmincon.ConstraintTolerance);


% Get fmincon options
nlpOptions = options.fmincon;

% load("matlab.mat" )
% constrFun_handle(decvar,constraintGrid)
% [~, grad_f] = objFun_handle(decvar)
% constraintGrid = subdvd_T;
% current_decision_vars = decvar;

% For storing iteration history (optional but good practice)
history.constraintGrid{iterationCount} = constraintGrid;
history.decisionVars{iterationCount} = current_decision_vars;

fprintf('\n--- Starting Interval Method Solver ---\n');

%% --- 2. Main Iteration Loop ---
while true
    fprintf('\n--- Iteration: %d ---\n', iterationCount);
    fprintf('Number of constraint grid points: %d\n', length(constraintGrid));

    % --- Call fmincon directly (replaces RPCDO_Solver.m) ---
    % Pass the dynamically changing constraintGrid to the constraint function
    % currentConstrFun = @(u) constrFun(u, constraintGrid);
    
    [decvar, obj_val, exitFlag, output, lambda, grad_f] = fmincon(...
        objFun_handle, current_decision_vars, [], [], [], [], lowerBounds, upperBounds, ...
        @(u)fminconConstrFun(u, constraintGrid), nlpOptions);
    
    history.cpuTime(iterationCount) = toc(startTime);
    
    if exitFlag > 0
        fprintf('NLP subproblem solved successfully.\n');
        current_decision_vars = decvar;
        history.isFeasible(iterationCount) = 1;
    else
        fprintf('Warning: NLP subproblem failed or was infeasible (Exit Flag: %d).\n', exitFlag);
        history.isFeasible(iterationCount) = 0;
    end
    
    % --- Check for Convergence ---
    t_midpoints = 0.5 * (constraintGrid(1:end-1) + constraintGrid(2:end));
    t_midpoints_cell = {t_midpoints};
    
    % Step B: Evaluate the discrete path constraints at these midpoints.
    [cieq_Dis, ~, gCieq_Dis, ~] = disConstrFun(decvar, t_midpoints_cell);

    % Step C: Identify active constraints based on their values at the midpoints.
    activePathIdx = find(cieq_Dis > -options.activeConstraintTol & cieq_Dis <= 0);
    gCieqActive = gCieq_Dis(:, activePathIdx);

    % Step D: Check the KKT stationarity condition.
    [isConverged, kktValue] = KKT_MinStatV3(decvar, grad_f, gCieqActive, problem, options);

    if isConverged
        fprintf('\nConvergence criteria met (KKT Norm: %.4e). Terminating.\n', kktValue);
        break;
    end


    % --- Interval Refinement Logic (Implemented in Main Solver) ---    
    if exitFlag < 0
        % STRATEGY 1: If NLP solver failed, refine all intervals.
        intervalMidpoints = 0.5 * (constraintGrid(1:end-1) + constraintGrid(2:end));
        constraintGrid = unique([constraintGrid, intervalMidpoints]);

    else
        % STRATEGY 2: If NLP solver succeeded, refine intervals where
        % the constraint value (at the midpoint) is close to being active.
        % This corresponds to the 'ActCNoD' logic in your VrFKKT_NewSbT.
        [cieq_Intvl, ~, ~, ~] = fminconConstrFun(decvar, constraintGrid);
        indicesToRefine = find(cieq_Intvl > -options.activeConstraintTol)'; % Ensure row vector
        % Call the dedicated function to perform the refinement.
        constraintGrid = refine_grid_Interval(constraintGrid, indicesToRefine, problem);
        

        if isempty(indicesToRefine)
             fprintf('KKT not met, but no active constraints found to refine. Terminating to avoid infinite loop.\n');
             break;
        end
        fprintf('Found %d intervals to refine based on active constraints.\n', length(indicesToRefine));
    end
    
    % --- Update State for Next Iteration ---
    iterationCount = iterationCount + 1;

    history.constraintGrid{iterationCount} = constraintGrid;
    history.decisionVars{iterationCount} = current_decision_vars;
    
    if iterationCount > 100 % Safety break
        fprintf('Warning: Maximum number of iterations (100) reached. Terminating.\n');
        break;
    end
end

%% --- 3. Package Results ---
fprintf('\n--- Algorithm Finished ---\n');
totalTime = toc(startTime);

solution.decisionVariables = decvar;
solution.objectiveValue = obj_val;
solution.exitFlag = exitFlag;
solution.lambda = lambda;
solution.gradObjective = grad_f;
solution.finalConstraintGrid = constraintGrid;

stats.totalIterations = iterationCount;
stats.cpuTime = totalTime;
stats.fminconOutput = output;
stats.history = history; % Include iteration history for analysis

fprintf('Total iterations: %d\n', stats.totalIterations);
fprintf('Total CPU time: %.2f seconds\n', stats.cpuTime);

end



function [c, ceq, grad_c, grad_ceq] = shifted_constraint_wrapper(u, cGrid, problem, shift_value)
    % Call the original upper-bound constraint function
    [c_original, ceq, grad_c, grad_ceq] = problem.functions.constraintsIntvl(u, cGrid, problem);
    
    % Apply the upward shift to the inequality constraint values.
    % The gradient remains unchanged because the derivative of a constant (shift_value) is zero.
    c = c_original + shift_value;
end