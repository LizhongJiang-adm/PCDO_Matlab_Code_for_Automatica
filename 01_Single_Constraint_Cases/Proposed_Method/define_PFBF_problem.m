function problem = define_PFBF_problem()
% = a=======================================================================
% Define the PFBF Optimal Control Problem
% =========================================================================
% Description:
% This function creates a struct containing all the necessary parameters,
% bounds, initial conditions, and function handles for the PFBF problem
% with constant control parameterization.
%
% Outputs:
%   problem - A struct with all problem-specific definitions.
% =========================================================================

%% Problem-Specific Function Handles
% These handles point to the dynamics, objective, and constraint functions
% for this specific problem (PFBF with constant control).
problem.functions.odeSensitivity   = @OdeSnstvSqnt_CstCtrl_PFBF;
problem.functions.odeUpperBound    = @OdeSnstvSqntUB_CstCtrl_PFBF;
problem.functions.objective        = @ObjtvSqnt_CstCtrl_PFBF;
problem.functions.constraintsDis   = @CnstrDisSqnt_CstCtrl_PFBF;
problem.functions.constraintsUB    = @CnstrUBSqnt_CstCtrl_PFBF;
problem.functions.pathDerivatives  = @GetDotPthCstr_PFBF;

%% Time and State Initialization
problem.time.T0 = 0;                  % Initial time
problem.time.TF = 40;                 % Final time
problem.state.x0 = [1, 0.2, 0.001, 250]; % Initial state vector

%% Control Discretization
problem.control.N = 40;               % Number of control intervals
problem.control.grid = linspace(problem.time.T0, problem.time.TF, problem.control.N + 1);

%% Decision Variable Definition
problem.numVars = problem.control.N * 1 * 1; % Total number of decision variables (u)
problem.lowerBounds =  0 * ones(1, problem.numVars); % Lower bounds for u
problem.upperBounds = 10 * ones(1, problem.numVars); % Upper bounds for u

% Initial Guess for Decision Variables
% For this problem, the specified starting point is a vector of all ones.
problem.initialGuess = ones(1, problem.numVars);

%% Path Constraints Definition
problem.pathConstraints.num = 1; % Number of path constraints
% Path constraint function handle: g(t,x,u) <= 0
problem.pathConstraints.handle{1} = @(t, x, u) x(:, 2) - 0.5;
problem.pathConstraints.constraintGrid = problem.control.grid;
end

