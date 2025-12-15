function problem = define_VDPO_problem_aBB()
% =========================================================================
% Define the Van der Pol Oscillator (VDPO) Problem for the aBB Method
% =========================================================================
% Description:
% This function creates a struct containing all the necessary parameters for
% the VDPO problem, specifically configured for the aBB Method solver.
% It defines only the control-related aspects of the decision variables.
%
% Outputs:
%   problem - A struct with all problem-specific definitions.
% =========================================================================

%% Problem-Specific Function Handles
problem.functions.odeSensitivity   = @OdeSnstvSqnt_CstCtrl_VDPO; % Assumed based on original script
problem.functions.objective        = @ObjtvSqnt_CstCtrl_VDPO;
problem.functions.constraints_aBB  = @Constraints_aBB_VDPO;
problem.functions.constraintsDis    = @CnstrDisSqnt_CstCtrl_VDPO; 

%% Time, State, and System Dimensions
problem.time.T0 = 0;                  % Initial time
problem.time.TF = 5;                  % Final time
problem.state.x0 = [0, 1, 0];         % Initial state vector [x1, x2, cost]
problem.system.dimensions = 3;        % Dimension of the ODE system

%% Control Discretization and Bounds
problem.control.N = 20;               % Number of control intervals
problem.control.grid = linspace(problem.time.T0, problem.time.TF, problem.control.N + 1);
problem.control.lowerBound = -0.3;    % Scalar lower bound for control u
problem.control.upperBound =  1.0;    % Scalar upper bound for control u

%% Decision Variable Definition (CONTROL VARIABLES ONLY)
problem.numVars = problem.control.N;
problem.lowerBounds_control = problem.control.lowerBound * ones(problem.numVars, 1);
problem.upperBounds_control = problem.control.upperBound * ones(problem.numVars, 1);

% Initial Guess for CONTROL variables ONLY
initial_u_value = 0.5 * (problem.control.lowerBound + problem.control.upperBound);
problem.initialGuess_control = initial_u_value * ones(problem.numVars, 1);

%% aBB Algorithm-Specific Parameters (Static)
problem.algorithm_params.myalpha = 10;    % Convexification parameter
problem.algorithm_params.tau = 1e-3;      % Smoothing parameter for phi_tau

%% Path Constraint Definition (for final plotting)
problem.pathConstraints.num = 1; % Number of path constraints
% Path constraint: -x1 - 0.4 <= 0
problem.pathConstraints.handle{1} = @(t, x, u) -x(:, 1) - 0.4;

end