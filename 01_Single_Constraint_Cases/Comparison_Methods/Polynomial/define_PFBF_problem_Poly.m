function problem = define_PFBF_problem_Poly()
% =========================================================================
% Define the PFBF Problem for the Polynomial Method
% =========================================================================
% Description:
% This function creates a struct containing all the necessary parameters for
% the PFBF problem, configured for the Polynomial Approximation solver.
% It follows the standardized format of the VDPO problem definition.
%
% Outputs:
%   problem - A struct with all problem-specific definitions.
% =========================================================================

%% Problem-Specific Function Handles
problem.functions.odeSensitivity   = @OdeSnstvSqnt_CstCtrl_PFBF;
problem.functions.objective        = @ObjtvSqnt_CstCtrl_PFBF;
problem.functions.constraintsPoly  = @Constraints_Poly_PFBF; % The new, standardized version
problem.functions.constraintsDis   = @CnstrDisSqnt_CstCtrl_PFBF;


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


%% Decision Variable Definition (Control Variables Only)
problem.numVars = problem.control.N;
problem.lowerBounds = problem.control.lowerBound * ones(1, problem.numVars);
problem.upperBounds = problem.control.upperBound * ones(1, problem.numVars);

% Initial Guess for CONTROL variables
% Based on original script str_pnt = ones(1,N);
problem.initialGuess = ones(1, problem.numVars);

%% Path Constraint Definition (for plotting and analysis)
problem.pathConstraints.num = 1;
% Path constraint: x2 - 0.5 <= 0
problem.pathConstraints.handle{1} = @(t, x, u) x(:, 2) - 0.5;
problem.functions.dotPathConstraint= @(t,dx) dx(:,2);

end