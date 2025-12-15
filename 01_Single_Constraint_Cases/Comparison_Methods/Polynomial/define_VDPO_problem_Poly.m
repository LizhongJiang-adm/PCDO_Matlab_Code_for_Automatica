function problem = define_VDPO_problem_Poly()
% =========================================================================
% Define the Van der Pol Oscillator (VDPO) Problem for the Polynomial Method
% =========================================================================
% Description:
% This function creates a struct containing all the necessary parameters for
% the VDPO problem, configured for the Polynomial Approximation solver.
%
% Outputs:
%   problem - A struct with all problem-specific definitions.
% =========================================================================

%% Problem-Specific Function Handles
problem.functions.odeSensitivity   = @OdeSnstvSqnt_CstCtrl_VDPO;
problem.functions.objective        = @ObjtvSqnt_CstCtrl_VDPO;
problem.functions.constraintsPoly  = @Constraints_Poly_VDPO; % Main constraint function
problem.functions.constraintsDis   = @CnstrDisSqnt_CstCtrl_VDPO;       % For KKT/convergence check


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

%% Decision Variable Definition (Control Variables Only)
problem.numVars = problem.control.N;
problem.lowerBounds = problem.control.lowerBound * ones(1, problem.numVars);
problem.upperBounds = problem.control.upperBound * ones(1, problem.numVars);

% Initial Guess for CONTROL variables
% Based on original script str_pnt = ones(1,N);
initial_u_value = 1;%0.5 * (problem.control.lowerBound + problem.control.upperBound);
problem.initialGuess = initial_u_value * ones(problem.numVars, 1);

%% Path Constraint Definition (for plotting and analysis)
problem.pathConstraints.num = 1;
% Path constraint: -x1 - 0.4 <= 0
problem.pathConstraints.handle{1} = @(t, x, u) -x(:, 1) - 0.4;
problem.functions.dotPathConstraint= @(t,dx) - dx(:,1);

end