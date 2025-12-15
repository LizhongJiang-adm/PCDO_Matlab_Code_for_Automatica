function [solution, stats] = solve_with_DisMethod(problem, options)
% =========================================================================
% Solve OCP using the Direct Discretization (Benchmark) Method
% =========================================================================
% Description:
% This function solves an optimal control problem by applying path constraints
% only at a fixed set of discrete time points. It performs a single,
% non-iterative call to the NLP solver (fmincon).
%
% Inputs:
%   problem - A struct containing the full problem definition.
%   options - A struct containing the NLP solver settings.
%
% Outputs:
%   solution - A struct containing the final optimization results.
%   stats    - A struct with statistics about the solution process.
% =========================================================================

%% --- 1. Initialization ---
startTime = tic;

% Extract parameters from input structs
lowerBounds = problem.lowerBounds;
upperBounds = problem.upperBounds;
initialGuess = problem.initialGuess;
controlGrid = problem.control.grid;

% Extract function handles
objFun_handle = @(u) problem.functions.objective(u, controlGrid, problem.state.x0);
disConstrFun = problem.functions.constraints;
nlpOptions = options.fmincon;

% --- Constraint Shifting Strategy ---
fminconConstrFun = @(u) shifted_constraint_wrapper_dis(u, problem, nlpOptions.ConstraintTolerance);

fprintf('\n--- Starting Direct Discretization Method Solver ---\n');
fprintf('Applying constraints at %d fixed points.\n', length(controlGrid));

%% --- 2. Solve the NLP Problem in a Single Run ---
[decvar, objVal, exitFlag, output, lambda, grad_f] = fmincon(...
    objFun_handle, initialGuess, [], [], [], [], lowerBounds, upperBounds, ...
    @(u) fminconConstrFun(u), nlpOptions);

if exitFlag > 0
    fprintf('NLP solver found a solution successfully.\n');
else
    fprintf('Warning: NLP solver did not find a solution (Exit Flag: %d).\n', exitFlag);
end

%% --- 3. Package Results ---
fprintf('\n--- Algorithm Finished ---\n');
totalTime = toc(startTime);

solution.decisionVariables = decvar;
solution.objectiveValue = objVal;
solution.exitFlag = exitFlag;
solution.lambda = lambda;
solution.gradObjective = grad_f;

% For this method, there's no iteration, but we can report the time.
stats.totalIterations = 1; % Only one NLP solve
stats.cpuTime = totalTime;
stats.fminconOutput = output;

fprintf('Total CPU time: %.2f seconds\n', stats.cpuTime);

end

%% --- Helper Function for Constraint Shifting ---
function [c, ceq, grad_c, grad_ceq] = shifted_constraint_wrapper_dis(u, problem, shift_value)
    % This wrapper calls the discrete constraint function.
    % The 't_dis' is passed in a cell array as expected by 'CnstrDis...' functions.
    [c_original, ceq, grad_c, grad_ceq] = problem.functions.constraints(u, problem.control.grid, problem.state.x0);
    
    % Apply the upward shift to the inequality constraints
    c = c_original + shift_value;
end