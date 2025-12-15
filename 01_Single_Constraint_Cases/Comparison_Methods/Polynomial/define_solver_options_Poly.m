function options = define_solver_options_Poly()
% =========================================================================
% Define Solver Options for the Polynomial Approximation Method
% =========================================================================
% Description:
% This function creates a struct containing all the necessary options and
% tolerances for the iterative polynomial approximation solver.
%
% Outputs:
%   options - A struct with all solver-specific settings.
% =========================================================================

%% Algorithm Convergence Tolerances
% Main tolerances for the outer loop convergence checks.
options.activeConstraintTol = 1e-3; % epsilon_act
options.stationarityTol = 1e-3;     % epsilon_sat

%% Algorithm Behavior Parameters
% Parameters that control the refinement strategy.
options.alpha = 1/4;   % Exponent for new point calculation
options.M_max = 10;    % Maximum number of subdivisions per interval
options.epsilon_app_factor = 0.5; % Factor for epsilon_app (epsilon_act/2)

%% NLP Solver Configuration (fmincon)
% Settings for the inner NLP subproblem solver.
nlpConstraintTol = options.activeConstraintTol * 1e-3;
% For consistency, we can set the optimality tolerance similarly.
nlpOptimalityTol = options.stationarityTol * 1e-3;


options.fmincon = optimoptions(@fmincon,...
    'Algorithm', 'active-set',...
    'SpecifyObjectiveGradient', true,...
    'SpecifyConstraintGradient', true,...
    'ConstraintTolerance', nlpConstraintTol,...
    'OptimalityTolerance', nlpOptimalityTol,...
    'MaxFunctionEvaluations', 10000,...
    'StepTolerance', 1e-15,...         % From original script
    'FunctionTolerance', 1e-15,...     % From original script
    'CheckGradients', false,...
    'Display', 'iter'); % Default display

end