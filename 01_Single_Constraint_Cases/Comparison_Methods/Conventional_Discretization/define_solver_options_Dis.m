function options = define_solver_options_Dis()
% =========================================================================
% Define Solver Options for the Direct Discretization Method
% =========================================================================
% Description:
% This function configures the NLP solver (fmincon) for the direct
% discretization benchmark method. This method solves the problem in a
% single NLP run without any iterative refinement.
%
% Outputs:
%   options - A struct containing the NLP solver settings.
% =========================================================================

%% NLP Solver Configuration (fmincon)
% This benchmark method directly calls fmincon once, so these are the
% only required options.

% Define the tolerances for the NLP solver.
% These can be set to standard high-precision values.
nlpConstraintTol = 1e-6;
nlpOptimalityTol = 1e-6;

% Configure the optimoptions object for fmincon.
options.fmincon = optimoptions(@fmincon,...
    'Algorithm', 'active-set',... % 'sqp' is often a good general-purpose choice
    'SpecifyObjectiveGradient', true,...
    'SpecifyConstraintGradient', true,...
    'ConstraintTolerance', nlpConstraintTol,...
    'OptimalityTolerance', nlpOptimalityTol,...
    'MaxFunctionEvaluations', 20000,...
    'Display', 'iter'); % Show iterative display for the single run

end