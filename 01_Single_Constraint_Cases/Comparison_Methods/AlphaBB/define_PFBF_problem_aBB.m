function problem = define_PFBF_problem_aBB()
% =========================================================================
% Define the PFBF Problem for the aBB Method
% =========================================================================
% Description:
% This function creates a struct containing all the necessary parameters for
% the PFBF problem, specifically configured for the aBB Method solver.
% It follows the standardized format of other aBB problem definitions.
%
% Outputs:
%   problem - A struct with all problem-specific definitions.
% =========================================================================

%% Problem-Specific Function Handles
problem.functions.odeSensitivity   = @OdeSnstvSqnt_CstCtrl_PFBF;
problem.functions.objective        = @ObjtvSqnt_CstCtrl_PFBF;
problem.functions.constraints_aBB  = @Constraints_aBB_PFBF;
problem.functions.constraintsDis    = @CnstrDisSqnt_CstCtrl_PFBF; 

%% Time, State, and System Dimensions
problem.time.T0 = 0;                  % Initial time
problem.time.TF = 40;                 % Final time
problem.state.x0 = [1, 0.2, 0.001, 250]; % Initial state vector
problem.system.dimensions = 4;        % Dimension of the ODE system

%% Control Discretization and Bounds
problem.control.N = 40;               % Number of control intervals
problem.control.grid = linspace(problem.time.T0, problem.time.TF, problem.control.N + 1);
problem.control.lowerBound = 0;       % Scalar lower bound for control u
problem.control.upperBound = 10;      % Scalar upper bound for control u

%% Decision Variable Definition (CONTROL VARIABLES ONLY)
problem.numVars = problem.control.N; % Maintained for potential backward compatibility
problem.numControlVars = problem.control.N;
problem.lowerBounds_control = problem.control.lowerBound * ones(problem.numControlVars, 1);
problem.upperBounds_control = problem.control.upperBound * ones(problem.numControlVars, 1);

% Initial Guess for CONTROL variables ONLY
initial_u_value = 0.5; % Based on original script u0 = 0.5*ones(N,1)
problem.initialGuess_control = initial_u_value * ones(problem.numControlVars, 1);

%% aBB Algorithm-Specific Parameters (Static)
problem.algorithm_params.myalpha = 0.5;   % Convexification parameter
problem.algorithm_params.tau = 2e-2;      % Smoothing parameter for phi_tau

%% Path Constraint Definition (for final plotting)
problem.pathConstraints.num = 1; % Number of path constraints
% Path constraint: x2 - 0.5 <= 0
problem.pathConstraints.handle{1} = @(t, x, u) x(:, 2) - 0.5;

end