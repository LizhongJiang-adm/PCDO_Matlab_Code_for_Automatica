function newConstraintGrid = UpdateCstrIntvlV3(currentConstraintGrid, indicesToDivide, fminconExitFlag)
% =========================================================================
% Update (Refine) the Constraint Grid (V3 - Logic Preserved)
% =========================================================================
% Description:
% This function refines the constraint grid by adding the midpoints of
% selected sub-intervals. It strictly preserves the logic of the original
% UpdateCstrIntvlV2 function.
%
% If the NLP solver (fmincon) failed, it refines ALL sub-intervals.
% Otherwise, it refines only the intervals specified by indicesToDivide.
%
% Inputs:
%   currentConstraintGrid - The vector of time points for the current grid.
%   indicesToDivide       - A vector of indices of the intervals to be refined.
%   fminconExitFlag       - The exit flag returned by the fmincon solver.
%
% Outputs:
%   newConstraintGrid - The updated, refined grid containing both old and
%                       new points.
% =========================================================================

% Check if fmincon failed to find a solution (e.g., infeasible)
if fminconExitFlag < 0
    % If NLP solver failed, override the indices to refine ALL intervals.
    fprintf('NLP solver failed (Flag=%d). Refining all %d intervals.\n', fminconExitFlag, length(currentConstraintGrid)-1);
    indicesToDivide = 1:(length(currentConstraintGrid) - 1);
end

% If there are no intervals to divide (either initially or after the check),
% simply return the original grid.
if isempty(indicesToDivide)
    newConstraintGrid = currentConstraintGrid;
    return;
end

% Calculate the midpoints of the intervals selected for division.
intervalMidpoints = 0.5 * (currentConstraintGrid(indicesToDivide) + currentConstraintGrid(indicesToDivide + 1));

% Combine the old grid with the new midpoints.
% 'unique' also sorts the grid, which is essential.
newConstraintGrid = unique([currentConstraintGrid, intervalMidpoints]);

end