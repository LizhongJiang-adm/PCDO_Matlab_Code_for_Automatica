function problem = define_JCB_problem_Poly()
% =========================================================================
% Define the JCB Problem for the Polynomial Method
% =========================================================================
% Description:
% This function creates a struct containing all the necessary parameters for
% the JCB problem, configured for the Polynomial Approximation solver.
% It follows the standardized format of other Poly problem definitions.
%
% Outputs:
%   problem - A struct with all problem-specific definitions.
% =========================================================================

%% Problem-Specific Function Handles
problem.functions.odeSensitivity   = @OdeSnstvSqnt_CstCtrl_JCB;
problem.functions.objective        = @ObjtvSqnt_CstCtrl_JCB;
problem.functions.constraintsPoly  = @Constraints_Poly_JCB; % The new, standardized version
problem.functions.constraintsDis   = @CnstrDisSqnt_CstCtrl_JCB;
problem.functions.dotPathConstraint= @(t,dx) dx(:,2) - 16*(t - 0.5);

%% Time, State, and System Dimensions
problem.time.T0 = 0;                  % Initial time
problem.time.TF = 1;                  % Final time
problem.state.x0 = [0, -1];        % Initial state vector [x1, x2, cost]
problem.system.dimensions = 2;        % Dimension of the ODE system (from original sys_dms=3)

%% Control Discretization and Bounds
problem.control.N = 20;               % Number of control intervals
problem.control.grid = linspace(problem.time.T0, problem.time.TF, problem.control.N + 1);
problem.control.lowerBound = -20;     % Scalar lower bound for control u
problem.control.upperBound =  20;     % Scalar upper bound for control u
problem.control.class = "Constant";   % For use by PathCnstrProfileV2

%% Decision Variable Definition (Control Variables Only)
problem.numVars = problem.control.N;
problem.lowerBounds = problem.control.lowerBound * ones(1, problem.numVars);
problem.upperBounds = problem.control.upperBound * ones(1, problem.numVars);

% Initial Guess for CONTROL variables
% Based on original script str_pnt = ones(1,N);
problem.initialGuess = ones(1, problem.numVars);

%% Path Constraint Definition (for plotting and analysis)
problem.pathConstraints.num = 1;
% Path constraint: x2 - 8*(t - 0.5)^2 + 0.5 <= 0
problem.pathConstraints.handle{1} = @(t, x, u) x(:, 2) - 8*(t - 0.5).^2 + 0.5;

end