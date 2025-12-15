function [solution, stats] = solve_with_UBMethod(problem, options)
% =========================================================================
% Solve Optimal Control Problem using the Iterative Upper-Bound Method
% =========================================================================
% Description:
% This function implements the iterative algorithm to handle path constraints
% by refining an upper-bound approximation. It calls a standard NLP
% solver (fmincon) in each iteration.
%
% Inputs:
%   problem - A struct containing the full problem definition (dynamics,
%             bounds, constraints, etc.).
%   options - A struct containing solver settings and tolerances.
%
% Outputs:
%   solution - A struct containing the final optimization results.
%   stats    - A struct with statistics about the solution process.
% =========================================================================

%% --- 1. Initialization ---

% %% REFAC: Start timer using tic/toc for better practice.
startTime = tic;
iterationCount = 1;

% %% REFAC: Get initial guess and bounds from the 'problem' struct.
initialGuess =  problem.initialGuess;
lowerBounds = problem.lowerBounds;
upperBounds = problem.upperBounds;

% %% REFAC: Initialize the constraint grid. This is a state of the solver.
% The initial grid for constructing upper-bound constraints.
constraintGrid = problem.pathConstraints.constraintGrid;

% %% REFAC: Define function handles using the 'problem' struct.
% This makes the solver independent of function names.
objFun = @(u) problem.functions.objective(u, problem.control.grid, problem.state.x0);
ubConstrFun = @(u, cGrid) problem.functions.constraintsUB(u, problem.control.grid, problem.state.x0, cGrid);
% To ensure the final solution strictly satisfies g(u) <= 0 rather than g(u) <= Tol,
% we shift the constraints passed to fmincon by its own tolerance.
disConstrFun = @(u, tMax) problem.functions.constraintsDis(u, problem.control.grid, problem.state.x0, tMax);

% %% REFAC: Make a local copy of the NLP options for fmincon.
nlpOptions = options.fmincon;

% Use a persistent variable for the starting point to enable warm starts.
currentDecisionVars = initialGuess;
fminconConstrFun = @(u, cGrid) shifted_constraint_wrapper(u, cGrid, problem, options.fmincon.ConstraintTolerance);


fprintf('\n--- Starting Upper-Bound Method ---\n');

%% --- 2. Main Iteration Loop ---
while true
    fprintf('\n--- Iteration: %d ---\n', iterationCount);
    fprintf('Number of constraint grid points: %d\n', length(constraintGrid));

    % Solve the NLP subproblem with the current upper-bound constraint grid
    [decisionVars, objValue, exitFlag, output, lambda, gradObj] = fmincon(...
        objFun, currentDecisionVars, [], [], [], [], lowerBounds, upperBounds, ...
        @(u) fminconConstrFun(u, constraintGrid), nlpOptions);

    history.CnstrGrid{iterationCount} = constraintGrid;
    history.iterationPoint{iterationCount} = decisionVars;
    history.objVal{iterationCount} = objValue;

    if exitFlag > 0 % fmincon found a solution
        % --- Check for convergence ---

        % Find the true maximum of the path constraint for the current solution
        [gMaxIntvl, tMaxIntvl] = MaxPathCstr_SubIntvlV3(decisionVars, constraintGrid, problem);
        tMaxIntvlCell = {tMaxIntvl};

        % Evaluate the discrete path constraints at these true maximum points
        [cieq_Dis, ~, gCieq_Dis, ~] = disConstrFun(decisionVars, tMaxIntvlCell);

        % %% REFAC: Use tolerances from the 'options' struct.
        activePathIdx = find(cieq_Dis > -options.activeConstraintTol & cieq_Dis <= 0);
        activePathTimes = tMaxIntvl(activePathIdx);
        gCieqActive = gCieq_Dis(:, activePathIdx);

        % Check KKT stationarity condition
        [isConverged, kktValue] = KKT_MinStatV3(decisionVars, gradObj, gCieqActive, problem, options); % Assuming KKT function can take tolerance
        if isConverged
            fprintf('\nConvergence criteria met. KKT conditions satisfied.\n');
            break;
        end

        % --- Identify intervals to refine ---
        [cieq_UB, ~, gCieq_UB, ~] = fminconConstrFun(decisionVars, constraintGrid);
        activeUBIdx = find(cieq_UB > -nlpOptions.ConstraintTolerance); % Active upper-bound constraints

        % Condition 1: Poor approximation of the max value
        % %% REFAC: Use tolerance from 'options'
        idx_value_mismatch = find(cieq_UB - gMaxIntvl > 0.5 * options.activeConstraintTol);
        
        % Condition 2: Poor approximation of the gradient
        % %% REFAC: Use tolerance from 'options'
        grad_norm_diff = vecnorm(gCieq_Dis - gCieq_UB);
        idx_grad_mismatch = find(grad_norm_diff > 0.5 * options.stationarityTol / max(1, numel(activeUBIdx)));

        % Combine indices of intervals that need refinement
        intervalsToRefine_Val = intersect(idx_value_mismatch, activeUBIdx);
        intervalsToRefine_Grad = intersect(idx_grad_mismatch, activeUBIdx);
        intervalsToRefine = unique([intervalsToRefine_Val, intervalsToRefine_Grad]);
        intervalsToRefine = find(cieq_UB > -2 * options.fmincon.ConstraintTolerance);

        % %% REFAC: Use current solution as the starting point for the next iteration (warm start).
        currentDecisionVars = decisionVars;

    else % fmincon failed
        fprintf('Warning: fmincon failed to find a solution (Exit Flag: %d). Terminating.\n', exitFlag);
        % Decide how to handle this. For now, we terminate.
        break;
    end



    if isempty(intervalsToRefine)
        fprintf('No intervals to refine, but KKT not met. Check tolerances or problem formulation. Terminating.\n');
        break;
    end
    
    % --- Refine the constraint grid ---
    fprintf('Refining %d interval(s).\n', numel(intervalsToRefine));

    newPoints = UpdateCstrIntvlV3(constraintGrid, intervalsToRefine, exitFlag);
    constraintGrid = unique(newPoints);

    if iterationCount == 1
        plot_pathconstraint(decisionVars, problem);
        hold on
        title("Path Constraint of PFBF Problem at First Iteration")
        % yline(cieq_UB,'g',LineWidth=1.2)
        hold off
        exportgraphics(gcf, "OneCnstrGrid_PathConstraintPFBF_Iteration1.pdf", 'Resolution', 600, 'ContentType', 'vector');
    end
    

    iterationCount = iterationCount + 1;

    

    if iterationCount > 100 % Safety break
        fprintf('Warning: Maximum number of iterations (100) reached. Terminating.\n');
        break;
    end
end

%% --- 3. Package Results ---
fprintf('\n--- Algorithm Finished ---\n');
totalTime = toc(startTime);

% %% REFAC: Package results into nice structs for output.
solution.decisionVariables = decisionVars;
solution.objectiveValue = objValue;
solution.exitFlag = exitFlag;
solution.lambda = lambda;
solution.gradObjective = gradObj;
solution.finalConstraintGrid = constraintGrid;

stats.totalIterations = iterationCount;
stats.cpuTime = totalTime;
stats.fminconOutput = output;
stats.history = history;

fprintf('Total iterations: %d\n', stats.totalIterations);
fprintf('Total CPU time: %.2f seconds\n', stats.cpuTime);

end


function [c, ceq, grad_c, grad_ceq] = shifted_constraint_wrapper(u, cGrid, problem, shift_value)
    % Call the original upper-bound constraint function
    [c_original, ceq, grad_c, grad_ceq] = problem.functions.constraintsUB(u, problem.control.grid, problem.state.x0, cGrid);
    
    % Apply the upward shift to the inequality constraint values.
    % The gradient remains unchanged because the derivative of a constant (shift_value) is zero.
    c = c_original + shift_value;
end