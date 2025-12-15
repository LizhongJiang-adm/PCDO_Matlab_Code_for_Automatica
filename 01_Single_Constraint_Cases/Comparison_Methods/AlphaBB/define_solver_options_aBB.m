function options = define_solver_options_aBB()
% =========================================================================
% Define Solver Options for the aBB Method
% =========================================================================
% Description:
% This function creates a struct containing all the necessary options and
% tolerances for the iterative aBB solver algorithm. It specifically
% configures the underlying NLP solver (fmincon) with settings found in
% the original aBB implementation.
%
% Outputs:
%   options - A struct with all solver-specific settings.
% =========================================================================

%% Algorithm Convergence Tolerances
% These tolerances are for the main, outer loop of the aBB method.

% Tolerance for determining if a constraint is 'active' for refinement purposes.
% Corresponds to epsilon_act in the original code.
options.activeConstraintTol = 1e-3;

% Tolerance for the KKT stationarity condition check.
% Corresponds to epsilon_sat in the original code.
options.stationarityTol = 1e-3;

%% NLP Solver Configuration (fmincon)
% These settings are for the inner NLP subproblem solver.

% Set the tolerances for fmincon based on the main algorithm's tolerances,
% as specified in the original 'subproblem_aBB.m'.
nlpConstraintTol = 1e-3 * options.activeConstraintTol;
nlpOptimalityTol = 1e-3 * options.stationarityTol;

% Configure the optimoptions object for fmincon.
options.fmincon = optimoptions(@fmincon,...
    'Algorithm', 'sqp',...
    'SpecifyObjectiveGradient', true, ...
    'SpecifyConstraintGradient', true,...
    'ConstraintTolerance', nlpConstraintTol,...
    'OptimalityTolerance', nlpOptimalityTol,...
    'MaxFunctionEvaluations', 10000,... % From original script
    'Display', 'iter'); % Default display, can be overridden for timing.

end