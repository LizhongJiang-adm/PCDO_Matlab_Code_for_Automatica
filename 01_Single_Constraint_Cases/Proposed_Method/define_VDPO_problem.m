function problem = define_VDPO_problem()
% =========================================================================
% Define the Van der Pol Oscillator (VDPO) Optimal Control Problem
% =========================================================================
% Description:
% This function creates a struct containing all the necessary parameters,
% bounds, initial conditions, and function handles for the VDPO problem
% with constant control parameterization.
%
% Outputs:
%   problem - A struct with all problem-specific definitions.
% =========================================================================

%% Problem-Specific Function Handles
% These handles point to the dynamics, objective, and constraint functions
% for this specific problem (VDPO with constant control).
problem.functions.odeSensitivity   = @OdeSnstvSqnt_CstCtrl_VDPO;
problem.functions.odeUpperBound    = @OdeSnstvSqntUB_CstCtrl_VDPO;
problem.functions.objective        = @ObjtvSqnt_CstCtrl_VDPO;
problem.functions.constraintsDis   = @CnstrDisSqnt_CstCtrl_VDPO;
problem.functions.constraintsUB    = @CnstrUBSqnt_CstCtrl_VDPO;
problem.functions.pathDerivatives  = @GetDotPthCstr_VDPO;

%% Time and State Initialization
problem.time.T0 = 0;              % Initial time
problem.time.TF = 5;              % Final time
problem.state.x0 = [0, 1, 0];     % Initial state vector [x1, x2, cost]

%% Control Discretization
problem.control.N = 20;           % Number of control intervals
problem.control.grid = linspace(problem.time.T0, problem.time.TF, problem.control.N + 1);

%% Decision Variable Definition
problem.numVars = problem.control.N * 1 * 1; % Total number of decision variables (u)
problem.lowerBounds = -0.3 * ones(1, problem.numVars); % Lower bounds for u
problem.upperBounds =  1.0 * ones(1, problem.numVars); % Upper bounds for u

%% Path Constraints Definition
problem.pathConstraints.num = 1; % Number of path constraints
% Path constraint function handle: g(t,x,u) <= 0
problem.pathConstraints.handle{1} = @(t, x, u) -x(:, 1) - 2/5;
problem.pathConstraints.constraintGrid = linspace(problem.time.T0, problem.time.TF, problem.control.N + 1);


% Initial Guess for Decision Variables 
problem.initialGuess = 0.5 * (problem.lowerBounds + problem.upperBounds);

end