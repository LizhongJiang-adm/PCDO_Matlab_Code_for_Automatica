function problem = define_JCB_problem_Interval()
% =========================================================================
% Define the JCB Problem for the Interval Method
% =========================================================================
% Description:
% This function creates a struct containing all the necessary parameters
% for the JCB problem, configured for the Interval Method solver.
%
% Outputs:
%   problem - A struct with all problem-specific definitions.
% =========================================================================

%% Problem-Specific Function Handles
problem.functions.odeSensitivity    = @OdeSnstvSqnt_CstCtrl_JCB;
problem.functions.objective         = @ObjtvSqnt_CstCtrl_JCB;
problem.functions.constraintsIntvl  = @Constraints_Intvl_JCB;
problem.functions.constraintsDis    = @CnstrDisSqnt_CstCtrl_JCB; 

%% Time, State, and System Dimensions
problem.time.T0 = 0;                  % Initial time
problem.time.TF = 1;                  % Final time
problem.state.x0 = [0, -1];           % Initial state vector [x1, x2]
problem.system.dimensions = 2;        % Dimension of the ODE system

%% Control Discretization and Bounds
problem.control.N = 20;               % Number of control intervals
problem.control.grid = linspace(problem.time.T0, problem.time.TF, problem.control.N + 1);
problem.control.lowerBound = -15;     % Scalar lower bound for control u
problem.control.upperBound =  15;     % Scalar upper bound for control u

%% Decision Variable Definition
problem.numVars = problem.control.N;
problem.lowerBounds = problem.control.lowerBound * ones(problem.numVars, 1);
problem.upperBounds = problem.control.upperBound * ones(problem.numVars, 1);

% Initial Guess for Decision Variables
initial_u_value = 0.5 * (problem.control.lowerBound + problem.control.upperBound);
problem.initialGuess = initial_u_value * ones(problem.numVars, 1);

%% Algorithm-Specific Parameters for this Problem
% NOTE: These values might need to be tuned specifically for the JCB problem.
problem.algorithm_params.UBdot2g = 14; % Placeholder: Upper bound on g_ddot
problem.algorithm_params.useAcceleratedDivision = true; % Flag for accelerated division

%% Path Constraint Definition (for plotting)
problem.pathConstraints.num = 1; % Number of path constraints
% Path constraint: x2 - 8*(t - 0.5)^2 + 0.5 <= 0
problem.pathConstraints.handle{1} = @(t, x, u) x(:, 2) - 8.*(t - 0.5).^2 + 0.5;

end