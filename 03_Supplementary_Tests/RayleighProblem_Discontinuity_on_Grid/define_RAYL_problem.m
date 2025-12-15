function problem = define_RAYL_problem()
% = a=======================================================================
% Define the RAYL Optimal Control Problem
% =========================================================================
% Description:
% This function creates a struct containing all the necessary parameters,
% bounds, initial conditions, and function handles for the RAYL problem
% with constant control parameterization.
%
% Outputs:
%   problem - A struct with all problem-specific definitions.
% =========================================================================

%% Problem-Specific Function Handles
% These handles point to the dynamics, objective, and constraint functions
% for this specific problem (RAYL with constant control).
problem.functions.odeSensitivity   = @OdeSnstvSqnt_CstCtrl_RAYL;
problem.functions.odeUpperBound    = @OdeSnstvSqntUB_CstCtrl_RAYL;
problem.functions.objective        = @ObjtvSqnt_CstCtrl_RAYL;
problem.functions.constraintsDis   = @CnstrDisSqnt_CstCtrl_RAYL;
problem.functions.constraintsUB    = @CnstrUBSqnt_CstCtrl_RAYL;
problem.functions.pathDerivatives  = @GetDotPthCstr_RAYL;

%% Time and State Initialization
problem.time.T0 = 0;                  % Initial time
problem.time.TF = 4.5;                 % Final time
problem.state.x0 = [-5, -5, 0]; % Initial state vector

%% Control Discretization
problem.control.N = 20;               % Number of control intervals
problem.control.grid = linspace(problem.time.T0, problem.time.TF, problem.control.N + 1);

%% Decision Variable Definition
problem.numVars = problem.control.N * 1 * 1; % Total number of decision variables (u)
problem.lowerBounds =  -4 * ones(1, problem.numVars); % Lower bounds for u
problem.upperBounds = 4 * ones(1, problem.numVars); % Upper bounds for u

% Initial Guess for Decision Variables
% For this problem, the specified starting point is a vector of all ones.
problem.initialGuess = 0*ones(1, problem.numVars);

%% Path Constraints Definition
problem.pathConstraints.num = 1; % Number of path constraints
% Path constraint function handle: g(t,x,u) <= 0
problem.pathConstraints.handle{1} = @(t, x, u) u(:,1) + x(:,1)./6;
problem.pathConstraints.constraintGrid = problem.control.grid;
end

