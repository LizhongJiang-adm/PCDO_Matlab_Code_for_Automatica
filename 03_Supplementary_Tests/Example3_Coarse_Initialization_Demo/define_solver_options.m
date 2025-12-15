function options = define_solver_options()
% =========================================================================
% Define Solver Options for the Upper-Bound Method
% =========================================================================
% Description:
% This function creates a struct containing all the necessary options and
% tolerances for the iterative upper-bound solver algorithm. It also
% configures the underlying NLP solver (fmincon).
%
% Outputs:
%   options - A struct with all solver-specific settings.
% =========================================================================

%% Algorithm Tolerances
% Tolerance for determining if a path constraint is active at its maximum point.
% Corresponds to epsilon_act in the paper.
options.activeConstraintTol = 1e-3;

% Tolerance for the KKT stationarity condition.
% Corresponds to epsilon_sat in the paper.
options.stationarityTol = 1e-3;

%% Upper-Bound Method Strategy
% Grid for constructing the initial set of upper-bound constraints.
% Let's keep it tied to the control grid for the start.
% Note: This grid will be refined iteratively by the solver.
% We will initialize it inside the solver function itself.

% Heuristic for refining the constraint grid. If true, uses derivative
% information to add new grid points. If false, only subdivides intervals.
options.useHeuristicDivision = true;

%% NLP Solver Configuration (fmincon)
% Set tolerances for the NLP solver to be tighter than the main algorithm's
% tolerances to ensure the subproblems are solved accurately.
nlpConstraintTol = options.activeConstraintTol * 1e-3;
nlpOptimalityTol = options.stationarityTol * 1e-3;

options.fmincon = optimoptions(@fmincon,...
    'Algorithm', 'active-set',...
    'SpecifyObjectiveGradient', true,...
    'SpecifyConstraintGradient', true,...
    'ConstraintTolerance', nlpConstraintTol,...
    'OptimalityTolerance', nlpOptimalityTol,...
    'MaxFunctionEvaluations', 10000,...
    'CheckGradients', false,...
    'Display', 'none'); % 'iter' provides useful feedback during runtime

end