function plot_horizontal_lines(intervals, values, varargin)
% =========================================================================
% Plot Horizontal Lines over Specified Intervals
% =========================================================================
% Description:
% This function plots a series of horizontal lines. For each interval
% defined by two consecutive points in the 'intervals' vector, it draws a
% horizontal line segment at the height specified by the corresponding
% entry in the 'values' vector.
%
% This is an alternative to 'stairs' when the vertical connecting lines
% are not desired.
%
% Inputs:
%   intervals - A vector of interval endpoints of size (1, M+1) or (M+1, 1).
%               e.g., [t0, t1, t2, ..., tM].
%   values    - A vector of values for each interval of size (1, M) or (M, 1).
%               value(i) is the height for the interval [intervals(i), intervals(i+1)].
%   varargin  - Optional additional arguments to be passed to the plot
%               function, e.g., 'Color', 'r', 'LineWidth', 2.
%
% Example:
%   t = 0:0.2:1;      % intervals = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
%   y = rand(1, 5);   % values for each of the 5 intervals
%   figure;
%   plot_horizontal_lines(t, y, 'Color', 'blue', 'LineWidth', 2, 'DisplayName', 'My Lines');
%   legend;
% =========================================================================

% --- Input Validation ---
if (length(intervals) - 1) ~= length(values)
    error('Input size mismatch: length(intervals) must be length(values) + 1.');
end

% Ensure inputs are row vectors for consistent indexing
intervals = intervals(:)';
values = values(:)';

% --- Plotting Loop ---
hold_state = ishold; % Check if hold is already on
hold on;

% The main plot function 'plot' takes pairs of X and Y coordinates.
% To draw a line from (x1, y) to (x2, y), we need X = [x1, x2] and Y = [y, y].
for i = 1:length(values)
    % Define the X and Y coordinates for the i-th horizontal line segment
    x_coords = [intervals(i), intervals(i+1)];
    y_coords = [values(i), values(i)];
    
    % Plot the segment
    plot(x_coords, y_coords, varargin{:});
    
    % --- Handle DisplayName for legend ---
    % The DisplayName is passed in varargin. To avoid creating multiple
    % legend entries, we turn off HandleVisibility for subsequent plots.
    if i==1
        % Find the DisplayName and remove it for the next iterations
        idx = find(strcmpi(varargin, "DisplayName"));
        varargin(idx:idx+1) = []; % Remove both 'DisplayName' and its value
        
        % Also turn off HandleVisibility for future lines in this loop
        varargin = [varargin, {'HandleVisibility', 'off'}];
    end
end

% Restore the original hold state if it was off
if ~hold_state
    hold off;
end

end