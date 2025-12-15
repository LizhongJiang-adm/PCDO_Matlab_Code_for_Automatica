function problem = define_JCB_problem_aBB()
% =========================================================================
% Define the JCB Problem for the aBB Method
% =========================================================================
% Description:
% This function creates a struct containing all the necessary parameters for
% the JCB problem, specifically configured for the aBB Method solver.
% It follows the standardized format of the VDPO problem definition.
%
% Outputs:
%   problem - A struct with all problem-specific definitions.
% =========================================================================

%% Problem-Specific Function Handles
problem.functions.odeSensitivity   = @OdeSnstvSqnt_CstCtrl_JCB;
problem.functions.objective        = @ObjtvSqnt_CstCtrl_JCB;
problem.functions.constraints_aBB  = @Constraints_aBB_JCB;
problem.functions.constraintsDis    = @CnstrDisSqnt_CstCtrl_JCB; 

%% Time, State, and System Dimensions
problem.time.T0 = 0;                  % Initial time
problem.time.TF = 1;                  % Final time
problem.state.x0 = [0, -1];        % Initial state vector 
problem.system.dimensions = 2;        % Dimension of the ODE system

%% Control Discretization and Bounds
problem.control.N = 20;               % Number of control intervals
problem.control.grid = linspace(problem.time.T0, problem.time.TF, problem.control.N + 1);
problem.control.lowerBound = -15;    % Scalar lower bound for control u
problem.control.upperBound =  15;    % Scalar upper bound for control u

%% Decision Variable Definition (CONTROL VARIABLES ONLY)
problem.numVars = problem.control.N; % This name is kept for backward compatibility if needed
problem.numControlVars = problem.control.N;
problem.lowerBounds_control = problem.control.lowerBound * ones(problem.numControlVars, 1);
problem.upperBounds_control = problem.control.upperBound * ones(problem.numControlVars, 1);

% Initial Guess for CONTROL variables ONLY
initial_u_value = 0; % Based on original script u0 = 0
problem.initialGuess_control = initial_u_value * ones(problem.numControlVars, 1);

%% aBB Algorithm-Specific Parameters (Static)
problem.algorithm_params.myalpha = 10;    % Convexification parameter
problem.algorithm_params.tau = 1e-2;      % Smoothing parameter for phi_tau

%% Path Constraint Definition (for final plotting)
problem.pathConstraints.num = 1; % Number of path constraints
% Path constraint: x2 - 8*(t - 0.5)^2 + 0.5 <= 0
problem.pathConstraints.handle{1} = @(t, x, u) x(:, 2) - 8.*(t - 0.5).^2 + 0.5;

end