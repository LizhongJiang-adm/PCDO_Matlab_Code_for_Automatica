


figure('Name', 'Algorithm Performance vs. Initial Grid Size', 'Position', [100, 100, 600, 450]);

% --- Define softer colors ---
soft_blue = [0.4, 0.6, 0.8];
soft_red = [0.84, 0.37, 0.33];

% --- Create the dual y-axis plot ---
yyaxis left; % Activate the left y-axis
plot(initial_interval_counts, cpu_times, '-o', 'LineWidth', 2, ...
    'Color', soft_blue, 'MarkerFaceColor', soft_blue, 'MarkerEdgeColor', soft_blue);
ylabel('\textbf{Total CPU Time (s)}', 'FontSize', 15,Interpreter='latex');
ax = gca;
ax.YColor = soft_blue; % Match axis color to the line color

tmp = cpu_times;
y_max = max(tmp);     y_min = min(tmp);
RangeU = y_max + 0.3*(y_max - y_min);
tmp = ylim;
RangeL = tmp(1);%y_min - 0.2*(y_max - y_min);
ylim([RangeL RangeU])


yyaxis right; % Activate the right y-axis
plot(initial_interval_counts, iteration_counts, '-s', 'LineWidth', 2, ...
    'Color', soft_red, 'MarkerFaceColor', soft_red, 'MarkerEdgeColor', soft_red,MarkerSize=8);
ylabel('\textbf{Number of Iterations}', 'FontSize', 15,Interpreter='latex');
ax = gca;
ax.YColor = soft_red; % Match axis color to the line color

tmp = iteration_counts;
y_max = max(tmp);     y_min = min(tmp);
RangeU = y_max + 0.2*(y_max - y_min);
RangeL = y_min-0.03;
ylim([RangeL RangeU])

% --- Common plot aesthetics ---
xlabel('\textbf{Number of Initial Subintervals}', 'FontSize', 17, Interpreter='latex');
xlim([1 1800])
legend('CPU Time', 'Number of Iterations', 'FontSize', 15, Interpreter='latex');
grid on;
box on;

set(gca, 'XScale', 'log');