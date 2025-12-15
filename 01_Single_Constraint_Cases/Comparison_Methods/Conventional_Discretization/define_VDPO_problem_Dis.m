function problem = define_VDPO_problem_Dis()
% =========================================================================
% Define the VDPO Problem for the Direct Discretization Method
% =========================================================================
% Description:
% This function creates a struct for the VDPO problem, configured for the
% direct discretization (benchmark) solver. Constraints are applied only at
% a fixed set of discrete points.
%
% Outputs:
%   problem - A struct with all problem-specific definitions.
% =========================================================================

%% Problem-Specific Function Handles
problem.functions.odeSensitivity   = @OdeSnstvSqnt_CstCtrl_VDPO;
problem.functions.objective        = @ObjtvSqnt_CstCtrl_VDPO;
% The constraint function is the simple discrete point constraint function.
problem.functions.constraints      = @Constraints_Dis_VDPO;

%% Time, State, and System Dimensions
problem.time.T0 = 0;                  % Initial time
problem.time.TF = 5;                  % Final time
problem.state.x0 = [0, 1, 0];         % Initial state vector [x1, x2, cost]
problem.system.dimensions = 3;        % Explicitly define state dimension

%% Control Discretization and Bounds
problem.control.N = 20;               % Number of control intervals
problem.control.grid = linspace(problem.time.T0, problem.time.TF, problem.control.N + 1);
problem.control.lowerBound = -0.3;    % Scalar lower bound for control u
problem.control.upperBound =  1.0;    % Scalar upper bound for control u

%% Decision Variable Definition
problem.numVars = problem.control.N;
problem.lowerBounds = problem.control.lowerBound * ones(1, problem.numVars);
problem.upperBounds = problem.control.upperBound * ones(1, problem.numVars);

% Initial Guess for CONTROL variables
initial_u_value = 1;
problem.initialGuess = initial_u_value * ones(problem.numVars, 1);

%% Path Constraint Definition (for plotting and analysis)
problem.pathConstraints.num = 1;
problem.pathConstraints.handle{1} = @(t, x, u) -x(:, 1) - 0.4;
% NOTE: dotPathConstraint is not needed for this method.

%% Discretization Grid for Constraints
% For this method, constraints are applied at the control grid points.
problem.constraintGrid = problem.control.grid;

end