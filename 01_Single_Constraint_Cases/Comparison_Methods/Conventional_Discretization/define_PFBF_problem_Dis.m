function problem = define_PFBF_problem_Dis()
% =========================================================================
% Define the PFBF Problem for the Direct Discretization Method
% =========================================================================
% ... (Description) ...

%% Problem-Specific Function Handles
problem.functions.odeSensitivity   = @OdeSnstvSqnt_CstCtrl_PFBF;
problem.functions.objective        = @ObjtvSqnt_CstCtrl_PFBF;
problem.functions.constraints      = @Constraints_Dis_PFBF;

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
problem.lowerBounds = problem.control.lowerBound * ones(1, problem.numVars);
problem.upperBounds = problem.control.upperBound * ones(1, problem.numVars);

% Initial Guess for CONTROL variables
problem.initialGuess = ones(1, problem.numVars);

%% Path Constraint Definition (for plotting and analysis)
problem.pathConstraints.num = 1;
problem.pathConstraints.handle{1} = @(t, x, u) x(:, 2) - 0.5;

%% Discretization Grid for Constraints
problem.constraintGrid = problem.control.grid;

end