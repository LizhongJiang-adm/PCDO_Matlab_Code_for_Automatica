function [solution, stats] = solve_with_PolyMethod(problem, options)
% =========================================================================
% Solve OCP using the Iterative Polynomial Approximation Method
% =========================================================================
% Description:
% This function implements an iterative algorithm for solving optimal
% control problems with path constraints using a polynomial approximation
% method. It calls a standard NLP solver (fmincon) in each iteration and
% refines a discrete time grid ('T_Dis') and adjusts residual parameters
% ('ResPmtr') until convergence criteria are met.
%
% Inputs:
%   problem - A struct containing the full problem definition.
%   options - A struct containing solver settings and tolerances.
%
% Outputs:
%   solution - A struct containing the final optimization results.
%   stats    - A struct with statistics about the solution process.
% =========================================================================

% =========================================================================
%                        IMPORTANT IMPLEMENTATION NOTE
%        Handling Discontinuities in Path Constraint Derivatives
% =========================================================================
% This optimal control problem uses piecewise-constant control inputs u(t).
% As a result, the time derivative of the path constraint, g_dot(t), is
% also a piecewise function and is DISCONTINUOUS at the control grid points
% (problem.control.grid).
%
% Any function that constructs a Hermite polynomial approximation of g(t)
% (e.g., Constraints_Poly_VDPO, GetPolyApproxError_PerInterval) MUST
% account for this. Specifically, when constructing a polynomial over an
% interval [t_a, t_b]:
%
% 1. The derivative at the left endpoint t_a MUST be the RIGHT-SIDED
%    derivative, g_dot(t_a+), calculated using the control u(t) that is
%    active for t >= t_a.
%
% 2. The derivative at the right endpoint t_b MUST be the LEFT-SIDED
%    derivative, g_dot(t_b-), calculated using the control u(t) that is
%    active for t < t_b.
%
% Failure to use the correct one-sided derivatives will result in large
% approximation errors for any polynomial interval that contains or ends
% at a control switching point.
% =========================================================================

%% --- 1. Initialization ---
startTime = tic;
iterationCount = 1;

% Extract parameters from input structs
lowerBounds = problem.lowerBounds;
upperBounds = problem.upperBounds;
current_decision_vars = problem.initialGuess;
controlGrid = problem.control.grid;

% Initialize the dynamic state of the solver
T_Dis = controlGrid; % Start with the control grid as the initial discrete grid
ResPmtr = zeros(1, length(T_Dis)-1); % Initialize residual parameters to zero

% Extract function handles
objFun_handle = @(u) problem.functions.objective(u, controlGrid, problem.state.x0);
polyConstrFun = problem.functions.constraintsPoly;
disConstrFun = @(u, tMax) problem.functions.constraintsDis(u, problem.control.grid, problem.state.x0, tMax);
nlpOptions = options.fmincon;

% --- Constraint Shifting Strategy ---
% To ensure strict feasibility g(u)<=0, we shift the constraints for fmincon.
fminconConstrFun = @(u, t_dis, res_pmtr) shifted_constraint_wrapper_poly(u, t_dis, res_pmtr, problem, nlpOptions.ConstraintTolerance);

% History tracking
history.T_Dis{iterationCount} = T_Dis;
history.ResPmtr{iterationCount} = ResPmtr;
history.decisionVars{iterationCount} = current_decision_vars;

fprintf('\n--- Starting Polynomial Approximation Method Solver ---\n');

%% --- 2. Main Iteration Loop ---
while true
    fprintf('\n--- Iteration: %d ---\n', iterationCount);
    fprintf('Number of discrete points (T_Dis): %d\n', length(T_Dis));


    % --- Call fmincon ---
    [decvar, objVal, exitFlag, output, lambda, grad_f] = fmincon(...
        objFun_handle, current_decision_vars, [], [], [], [], lowerBounds, upperBounds, ...
        @(u) fminconConstrFun(u, T_Dis, ResPmtr), nlpOptions);

    % load('matlab.mat')
    % load('matlab2.mat')
    % T_Dis = T_Dis_org;
    % decvar = DecVar_org;
    % ResPmtr = ResPmtr_org;    

    current_decision_vars = decvar;
    history.cpuTime(iterationCount) = toc(startTime);
    history.decisionVars{iterationCount+1} = current_decision_vars;

    % ================================================================
    % --- Convergence Check and Refinement (Polynomial Method Logic) ---
    % ================================================================

    % --- Setup for Checks ---
    % Create the full interval set for analysis, including t0 and tf
    poly_intervals = unique([problem.time.T0, T_Dis, problem.time.TF]);
    num_poly_intervals = length(poly_intervals) - 1;

    % Calculate the true maximum of the path constraint over the whole trajectory
    [g_max_total, ~] = PathCnstrProfileV3(decvar, problem); % Requires a refactored PathCnstrProfile
    % [g_max_total, t_at_max, state_and_g, time_vector, control_vector] = PathCnstrProfileV3(decvar, problem); % Requires a refactored PathCnstrProfile
    % plot(time_vector,state_and_g(:,end))

    % Calculate the max of the polynomial approximation in each interval
    TmaxPoly = GetPolyTmax_PerInterval_Poly(decvar, problem, T_Dis); % Placeholder for a refactored version

    % Re-evaluate the constraints to get their un-shifted values
    [cieq, ~, grad_cieq, ~] = fminconConstrFun(decvar, T_Dis, ResPmtr);
    num_T_Dis = length(T_Dis);
    CstrPoly = cieq(num_T_Dis+1:end);     
    GCstrPoly = grad_cieq(:,num_T_Dis+1:end);
    
    % Get Lagrange multipliers for the polynomial constraints
    LmbPoly = lambda.ineqnonlin(num_T_Dis+1:end);
    
    % Evaluate true constraint values at the polynomial maximum points
    [CieMaxItv, ~, GCieMaxItv, ~] = disConstrFun(decvar, {TmaxPoly});
    
    PolyApmVal = CstrPoly - ResPmtr; % Approximation value without residual
    
    % --- Check Termination Conditions ---
    IdxActPoly = find(CstrPoly > -options.activeConstraintTol * 0.5);
    
    val_error_check = all(abs(PolyApmVal(IdxActPoly) - CieMaxItv(IdxActPoly)) <= options.activeConstraintTol * 0.5);
    % grad_error_check = all(LmbPoly(IdxActPoly)' .* (GCstrPoly(:,IdxActPoly) - GCieMaxItv(:,IdxActPoly)) <= options.stationarityTol * 0.5 / max(1, length(IdxActPoly)));
    grad_error_check = all(vecnorm(LmbPoly(IdxActPoly)' .* (GCstrPoly(:,IdxActPoly) - GCieMaxItv(:,IdxActPoly))) <= options.stationarityTol * 0.5 / max(1, length(IdxActPoly)) );

    if g_max_total <= 0 && val_error_check && grad_error_check
        fprintf('\nConvergence criteria met. Terminating.\n');
        break;
    end
    
    % --- Identify Intervals to Refine (IdxErr) ---
    if g_max_total <= 0
        % Errors are due to poor approximation
        Tmp = find(abs(CstrPoly - CieMaxItv) > options.activeConstraintTol * 0.5);
        IdxErrVal = intersect(Tmp, IdxActPoly);
        grad_error_terms = LmbPoly' .* vecnorm(GCstrPoly - GCieMaxItv);
        IdxErrGra = find(grad_error_terms > options.stationarityTol * 0.5 / max(1, length(IdxActPoly)));
        IdxErr = unique([IdxErrVal, IdxErrGra]);
    else
        % Errors are due to constraint violation
        % [GmaxIntvl, ~] = MaxPathCstr_SubIntvlV3(decvar, poly_intervals, problem);
        [GmaxIntvl, TmaxIntvl] = GetPathConstraintMax_PerInterval(decvar, problem, poly_intervals);
        IdxVlt = find(GmaxIntvl > 0);
        IdxErr = IdxVlt;
    end
    
    % --- Update T_Dis and ResPmtr ---
    [PolyErr, ~] = GetPolyApproxError_PerInterval(decvar, problem, T_Dis); % Placeholder for refactored version

    ResPNew = cell(1, num_poly_intervals);
    current_res_pmtr_idx = 1;

    for i = 1:num_poly_intervals
        if ismember(i, IdxErr)
            % This interval needs refinement
            epsilon_app = options.activeConstraintTol * options.epsilon_app_factor;
            
            Wtd = poly_intervals(i+1) - poly_intervals(i);
            DeltT = Wtd * (epsilon_app / PolyErr(i))^options.alpha;
            mDvd = Wtd / DeltT;
            
            new_points_in_interval = GetNewPnt(poly_intervals(i:i+1), mDvd, options.M_max);
            T_Dis = sort([T_Dis, new_points_in_interval]);
            
            num_new_subintervals = length(new_points_in_interval) + 1;
            ResPNew{i} = epsilon_app * ones(1, num_new_subintervals);
        else
            % This interval is not refined, keep its ResPmtr
            ResPNew{i} = ResPmtr(current_res_pmtr_idx);
            current_res_pmtr_idx = current_res_pmtr_idx + 1;
        end
    end
    ResPmtr = [ResPNew{:}];

    % --- Update State for Next Iteration ---
    iterationCount = iterationCount + 1;
    history.T_Dis{iterationCount} = T_Dis;
    history.ResPmtr{iterationCount} = ResPmtr;
    
    if iterationCount > 100, break; end % Safety break
end

%% --- 3. Package Results ---
fprintf('\n--- Algorithm Finished ---\n');
solution.decisionVariables = decvar;
solution.objectiveValue = objVal;
solution.exitFlag = exitFlag;
solution.final_T_Dis = T_Dis;
solution.final_ResPmtr = ResPmtr;
% ... other stats ...
stats.totalIterations = iterationCount;
stats.cpuTime = toc(startTime);
stats.history = history;

fprintf('Total iterations: %d\n', stats.totalIterations);
fprintf('Total CPU time: %.2f seconds\n', stats.cpuTime);

end

%% --- Helper Functions ---
function [c, ceq, grad_c, grad_ceq] = shifted_constraint_wrapper_poly(u, t_dis, res_pmtr, problem, shift_value)
    [c_original, ceq, grad_c, grad_ceq] = problem.functions.constraintsPoly(u, problem, t_dis, res_pmtr);
    
    c = c_original + shift_value;
end

function NewPnt = GetNewPnt(End, mDvd, M_max)
    % (This is a direct copy of your original GetNewPnt function)
    if 2<mDvd && mDvd<=M_max, NumDvd = ceil(mDvd);
    elseif mDvd>M_max, NumDvd = M_max;
    else, NumDvd = 2;
    end
    LEnd = End(1); REnd = End(end);
    NewPnt = LEnd : (REnd - LEnd)/NumDvd : REnd;
    NewPnt([1, end]) = []; % Remove endpoints
end