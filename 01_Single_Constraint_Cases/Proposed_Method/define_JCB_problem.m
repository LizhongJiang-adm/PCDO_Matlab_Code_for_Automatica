function problem = define_JCB_problem()
% =========================================================================
% Define the JCB Optimal Control Problem
% =========================================================================
% Description:
% This function creates a struct containing all the necessary parameters,
% bounds, initial conditions, and function handles for the JCB problem
% with constant control parameterization.
%
% Outputs:
%   problem - A struct with all problem-specific definitions.
% =========================================================================

%% Problem-Specific Function Handles
% These handles point to the dynamics, objective, and constraint functions
% for this specific problem (JCB with constant control).
problem.functions.odeSensitivity   = @OdeSnstvSqnt_CstCtrl_JCB;
problem.functions.odeUpperBound    = @OdeSnstvSqntUB_CstCtrl_JCB;
problem.functions.objective        = @ObjtvSqnt_CstCtrl_JCB;
problem.functions.constraintsDis   = @CnstrDisSqnt_CstCtrl_JCB;
problem.functions.constraintsUB    = @CnstrUBSqnt_CstCtrl_JCB;
problem.functions.pathDerivatives  = @GetDotPthCstr_JCB;

%% Time and State Initialization
problem.time.T0 = 0;          % Initial time
problem.time.TF = 1;          % Final time
problem.state.x0 = [0, -1];   % Initial state vector [x1, x2]

%% Control Discretization
problem.control.N = 20;       % Number of control intervals
problem.control.grid = linspace(problem.time.T0, problem.time.TF, problem.control.N + 1);

%% Decision Variable Definition
problem.numVars = problem.control.N * 1 * 1; % Total number of decision variables (u)
problem.lowerBounds = -15 * ones(1, problem.numVars); % Lower bounds for u
problem.upperBounds =  15 * ones(1, problem.numVars); % Upper bounds for u

% Initial Guess for Decision Variables
% Using the midpoint of the bounds as a safe default.
problem.initialGuess = 0.5 * (problem.lowerBounds + problem.upperBounds);

%% Path Constraints Definition
problem.pathConstraints.num = 1; % Number of path constraints
% Path constraint function handle: g(t,x,u) <= 0
problem.pathConstraints.handle{1} = @(t, x, u) x(:, 2) - 8 .* (t - 0.5).^2 + 0.5;
problem.pathConstraints.constraintGrid = problem.control.grid;
end