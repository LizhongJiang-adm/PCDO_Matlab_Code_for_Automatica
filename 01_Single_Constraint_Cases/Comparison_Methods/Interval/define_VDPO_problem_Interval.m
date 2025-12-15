function problem = define_VDPO_problem_Interval()
% =========================================================================
% Define the Van der Pol Oscillator (VDPO) Problem for the Interval Method
% =========================================================================
% Description:
% This function creates a struct containing all the necessary parameters,
% bounds, initial conditions, and function handles for the VDPO problem,
% specifically configured for the Interval Method solver.
%
% Outputs:
%   problem - A struct with all problem-specific definitions.
% =========================================================================

%% Problem-Specific Function Handles
% These handles point to the dynamics, objective, and constraint functions
% for this specific problem (VDPO for Interval Method).
problem.functions.odeSensitivity    = @OdeSnstvSqnt_CstCtrl_VDPO;
problem.functions.objective         = @ObjtvSqnt_CstCtrl_VDPO;
problem.functions.constraintsIntvl  = @Constraints_Intvl_VDPO;
problem.functions.constraintsDis    = @CnstrDisSqnt_CstCtrl_VDPO; 

%% Time, State, and System Dimensions
problem.time.T0 = 0;                  % Initial time
problem.time.TF = 5;                  % Final time
problem.state.x0 = [0, 1, 0];         % Initial state vector [x1, x2, cost]
problem.system.dimensions = 3;        % Dimension of the ODE system (OdeDm)

%% Control Discretization and Bounds
problem.control.N = 20;               % Number of control intervals
problem.control.grid = linspace(problem.time.T0, problem.time.TF, problem.control.N + 1);
problem.control.lowerBound = -0.3;    % Scalar lower bound for control u
problem.control.upperBound =  1.0;    % Scalar upper bound for control u

%% Decision Variable Definition
% For constant control, one variable per interval.
problem.numVars = problem.control.N;
problem.lowerBounds = problem.control.lowerBound * ones(problem.numVars, 1);
problem.upperBounds = problem.control.upperBound * ones(problem.numVars, 1);

% Initial Guess for Decision Variables
% Set to the midpoint of the bounds.
initial_u_value = 0.5 * (problem.control.lowerBound + problem.control.upperBound);
problem.initialGuess = initial_u_value * ones(problem.numVars, 1);

% These parameters are specific to the interval algorithm's strategy for this problem.
problem.algorithm_params.UBdot2g = 40; % Upper bound on the second derivative of the path constraint.
problem.algorithm_params.useAcceleratedDivision = true; % Flag for accelerated division (AclDvd = 1).

% This section defines the actual path constraint g(t,x,u) <= 0,
% which is needed by the plotting function to generate the final profile.
problem.pathConstraints.num = 1; % Number of path constraints
% The path constraint for VDPO is -x1 - 2/5 <= 0
problem.pathConstraints.handle{1} = @(t, x, u) -x(:, 1) - 2/5;

end