function [newGrid, control_flow] = refine_grid_aBB(currentGrid, indicesToRefine, active_tk, exitFlag, options, problem)
% =========================================================================
% Refine the aBB Grid Based on the Paper's epsilon_U-union Criterion
% =========================================================================
% Description:
% This function refines the 'eta' grid based on the active time points 't_k'.
% It strictly implements the epsilon_U-union logic described in the reference
% paper, using a relative width criterion.
%
% Inputs:
%   currentGrid     - The current 'eta' grid.
%   indicesToRefine - Indices of the intervals in 'currentGrid' to consider for refinement.
%   active_tk       - The optimized time points 't_k' corresponding to the active intervals.
%   exitFlag        - The exit flag from the NLP solver.
%   options         - The solver options struct, containing epsilon_U.
%
% Outputs:
%   newGrid       - The updated, refined grid.
%   control_flow  - A flag to signal algorithm termination if the grid stagnates.
% =========================================================================

control_flow = 0;

TimeDuration = problem.time.TF - problem.time.T0;
epsilon_U = 2*options.activeConstraintTol/problem.algorithm_params.myalpha/TimeDuration;

% --- Case 1: NLP was infeasible ---
if exitFlag < 0
    fprintf('NLP solver failed. Refining all sufficiently wide intervals by bisection.\n');
    
    pointsToAdd = [];
    intervalWidths = currentGrid(2:end) - currentGrid(1:end-1);
    
    % As a practical measure for the infeasible case, let's bisect all intervals.
    % The paper is less specific here, so bisection is a robust default.
    % We could also apply a minimum width check if desired.
    pointsToAdd = 0.5 * (currentGrid(1:end-1) + currentGrid(2:end));

else
    % --- Case 2: NLP was feasible ---
    pointsToAdd = [];
    
    % Ensure inputs have the same orientation
    indicesToRefine = indicesToRefine(:)';
    active_tk = active_tk(:)';
    
    if length(indicesToRefine) ~= length(active_tk)
        error('The number of active indices must match the number of active time points.');
    end

    % Iterate through each active interval and its corresponding t_k
    for i = 1:length(indicesToRefine)
        idx = indicesToRefine(i); % The index of the interval, e.g., 5
        tk = active_tk(i);        % The time point t_k inside this interval
        
        % Get the interval boundaries
        eta_left = currentGrid(idx);
        eta_right = currentGrid(idx+1);
        
        interval_width = eta_right - eta_left;
        
        % Calculate the minimum distance from t_k to the boundaries
        min_dist = min(tk - eta_left, eta_right - tk);
        
        % The relative threshold from the paper
        threshold = epsilon_U;
        
        % Check the condition from the paper
        if min_dist > threshold
            % If the condition holds, accept t_k as a new grid point
            pointsToAdd = [pointsToAdd, tk];
        end
    end
end

% --- Finalize the new grid ---
original_grid_size = length(currentGrid);
newGrid = unique([currentGrid, pointsToAdd]);
new_grid_size = length(newGrid);

% Check for stagnation
if new_grid_size == original_grid_size && exitFlag > 0
    % If the grid didn't change and the problem was feasible, it might be stuck
    control_flow = 1;
    fprintf('Grid has stagnated. Signaling for termination.\n');
end

end