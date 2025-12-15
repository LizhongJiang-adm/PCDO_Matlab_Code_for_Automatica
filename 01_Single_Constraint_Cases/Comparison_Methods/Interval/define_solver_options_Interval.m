function options = define_solver_options_Interval()
% =========================================================================
% Define Solver Options for the Interval Method
% =========================================================================
% Description:
% This function creates a struct containing all the necessary options and
% tolerances for the iterative interval-based solver algorithm. It also
% configures the underlying NLP solver (fmincon).
%
% Outputs:
%   options - A struct with all solver-specific settings.
% =========================================================================

%% Algorithm Convergence Tolerances
% These tolerances are for the main, outer loop of the interval method.

% Tolerance for determining if a constraint is 'active'.
% Corresponds to epsilon_act in the original code.
options.activeConstraintTol = 1e-3;

% Tolerance for the KKT stationarity condition check.
% Corresponds to epsilon_sat in the original code.
options.stationarityTol = 1e-3;

%% NLP Solver Configuration (fmincon)
% These settings are for the inner NLP subproblem solver.

% Set the constraint tolerance for fmincon. It's often set tighter
% than the main algorithm's tolerance to ensure subproblems are solved
% accurately. Corresponds to ActTol in the original code.
nlpConstraintTol = options.activeConstraintTol * 1e-3;

% For consistency, we can set the optimality tolerance similarly.
nlpOptimalityTol = options.stationarityTol * 1e-3;

% Configure the optimoptions object for fmincon.
options.fmincon = optimoptions(@fmincon,...
    'Algorithm', 'active-set',...
    'SpecifyObjectiveGradient', true,...
    'SpecifyConstraintGradient', true,...
    'ConstraintTolerance', nlpConstraintTol,...
    'OptimalityTolerance', nlpOptimalityTol, ...
    'MaxFunctionEvaluations', 10000,...
    'Display', 'iter'); % Default display, can be overridden for timing.

end