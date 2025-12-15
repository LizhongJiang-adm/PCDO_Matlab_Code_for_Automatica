function newGrid = refine_grid_Interval(currentGrid, indicesToRefine, problem)
% =========================================================================
% Refine the Constraint Grid for the Interval Method
% =========================================================================
% Description:
% This function refines the time grid by adding new points to the intervals
% specified by 'indicesToRefine'. It supports both standard bisection and
% an accelerated division strategy.
%
% Inputs:
%   currentGrid     - The current vector of time points.
%   indicesToRefine - A vector of indices of the intervals to be refined.
%   problem         - The problem definition struct, containing algorithm
%                     parameters like UBdot2g and the division strategy.
%
% Outputs:
%   newGrid - The updated, refined grid.
% =========================================================================

% If there are no intervals to refine, return the grid unchanged.
if isempty(indicesToRefine)
    newGrid = currentGrid;
    return;
end

% Extract necessary algorithm parameters from the problem struct
useAcceleratedDivision = problem.algorithm_params.useAcceleratedDivision;
UBdot2g = problem.algorithm_params.UBdot2g;
epsilon_act = 1e-3; % This should ideally be from options, but let's keep it for now.

% Calculate properties of the current grid
intervalWidths = currentGrid(2:end) - currentGrid(1:end-1);
intervalMidpoints = 0.5 * (currentGrid(1:end-1) + currentGrid(2:end));

pointsToAdd = [];

if useAcceleratedDivision
    % --- Accelerated Division Logic ---
    minWidth = sqrt(epsilon_act / (0.125 * UBdot2g));
    
    for i = 1:length(indicesToRefine)
        idx = indicesToRefine(i); % The actual index of the interval
        widthRatio = intervalWidths(idx) / minWidth;
        
        a = currentGrid(idx);
        b = currentGrid(idx+1);

        if widthRatio > 4
            % Subdivide into multiple smaller intervals
            numNewIntervals = ceil(widthRatio);
            for j = 1:(numNewIntervals - 1)
                newPoint = a + j * (b - a) / numNewIntervals;
                pointsToAdd = [pointsToAdd, newPoint];
            end
        else
            % Just add the midpoint
            pointsToAdd = [pointsToAdd, intervalMidpoints(idx)];
        end
    end
else
    % --- Standard Bisection Logic ---
    pointsToAdd = intervalMidpoints(indicesToRefine);
end

% Combine the old grid with the new points and ensure it's sorted and unique
newGrid = unique([currentGrid, pointsToAdd]);

end