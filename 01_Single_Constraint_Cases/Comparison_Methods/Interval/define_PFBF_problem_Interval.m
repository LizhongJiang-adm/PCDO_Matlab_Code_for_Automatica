function problem = define_PFBF_problem_Interval()
% =========================================================================
% Define the PFBF Problem for the Interval Method
% =========================================================================
% Description:
% This function creates a struct containing all the necessary parameters
% for the PFBF problem, configured for the Interval Method solver.
% The structure of this file is standardized based on the VDPO definition.
%
% Outputs:
%   problem - A struct with all problem-specific definitions.
% =========================================================================

%% Problem-Specific Function Handles
problem.functions.odeSensitivity    = @OdeSnstvSqnt_CstCtrl_PFBF;
problem.functions.objective         = @ObjtvSqnt_CstCtrl_PFBF;
problem.functions.constraintsIntvl  = @Constraints_Intvl_PFBF;
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

%% Decision Variable Definition
problem.numVars = problem.control.N;
problem.lowerBounds = problem.control.lowerBound * ones(problem.numVars, 1);
problem.upperBounds = problem.control.upperBound * ones(problem.numVars, 1);

% Initial Guess for Decision Variables
% The original script for this problem used a vector of all ones.
problem.initialGuess = ones(problem.numVars, 1);

%% Algorithm-Specific Parameters for this Problem
% NOTE: These values might need to be tuned specifically for the PFBF problem.
problem.algorithm_params.UBdot2g = 0.5; % Placeholder: Upper bound on g_ddot
problem.algorithm_params.useAcceleratedDivision = true; % Flag for accelerated division

%% Path Constraint Definition (for plotting)
problem.pathConstraints.num = 1; % Number of path constraints
% Path constraint: x2 - 0.5 <= 0
problem.pathConstraints.handle{1} = @(t, x, u) x(:, 2) - 0.5;

end