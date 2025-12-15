function [isConverged, kktStationarityValue] = KKT_MinStatV3(decisionVars, gradObjective, gradActivePathCstr, problem, options)
% =========================================================================
% Check KKT Stationarity Condition for the Optimal Control Problem
% =========================================================================
% Description:
% This function checks if the Karush-Kuhn-Tucker (KKT) stationarity
% condition is satisfied within a given tolerance. It considers the
% gradients of the active path constraints and the active decision
% variable bounds.
%
% The condition is: grad(Objective) + sum(lambda_i * grad(Constraint_i)) = 0
% where lambda_i >= 0. This is formulated as a non-negative least squares
% problem to find the optimal Lagrange multipliers (lambda_i).
%
% Inputs:
%   decisionVars      - The current vector of decision variables (u).
%   gradObjective     - The gradient of the objective function w.r.t. u.
%   gradActivePathCstr- The gradients of the active path constraints w.r.t. u.
%   problem           - The problem definition struct (contains bounds and dimensions).
%   options           - The solver options struct (contains tolerances).
%
% Outputs:
%   isConverged         - A flag (1 if converged, 0 otherwise).
%   kktStationarityValue- The norm of the residual of the KKT stationarity equation.
% =========================================================================

%% --- 1. Get Parameters from Input Structs ---
% No more global variables!
numVars = problem.numVars;
lowerBounds = problem.lowerBounds_control;
upperBounds = problem.upperBounds_control;

activeTol = options.activeConstraintTol; % epsilon_act
stationarityTol = options.stationarityTol; % epsilon_sat

%% --- 2. Assemble the Matrix of Active Constraint Gradients ---

% Initialize with the gradients from the active path constraints
allActiveGradients = gradActivePathCstr;

% Check for active upper and lower bound constraints on the decision variables
for i = 1:numVars
    % Check for active upper bound
    if decisionVars(i) >= upperBounds(i) - activeTol
        % Create the gradient for the constraint: u(i) - UB(i) <= 0
        grad_u_UB = zeros(numVars, 1);
        grad_u_UB(i) = 1;
        allActiveGradients = [allActiveGradients, grad_u_UB];
    end
    
    % Check for active lower bound
    if decisionVars(i) <= lowerBounds(i) + activeTol
        % Create the gradient for the constraint: LB(i) - u(i) <= 0
        grad_u_LB = zeros(numVars, 1);
        grad_u_LB(i) = -1;
        allActiveGradients = [allActiveGradients, grad_u_LB];
    end
end

% Handle the case where there are no active constraints
if isempty(allActiveGradients)
    % If there are no active constraints, the KKT condition simplifies to
    % checking if the gradient of the objective is zero.
    kktStationarityValue = norm(gradObjective);
else
    %% --- 3. Solve for Lagrange Multipliers ---
    % We want to solve: gradObjective + allActiveGradients * lambda = 0
    % This is equivalent to: allActiveGradients * lambda = -gradObjective
    % Since lambda must be non-negative, we use lsqnonneg.
    [lambda, ~] = lsqnonneg(allActiveGradients, -gradObjective);

    % Calculate the residual of the KKT stationarity condition
    residual = allActiveGradients * lambda + gradObjective;
    kktStationarityValue = norm(residual);
end

%% --- 4. Display Results and Determine Convergence ---
fprintf('KKT stationarity check: Norm = %.4e (Tolerance = %.4e)\n', kktStationarityValue, stationarityTol);

if kktStationarityValue <= stationarityTol
    isConverged = 1;
else
    isConverged = 0;
end

end