% =========================================================================
% benchmark_grid_performance_Ex1_VDPO
% =========================================================================
% Description:
% This script analyzes the performance of the Upper-Bound Method by varying
% the number of initial constraint intervals. It measures the total CPU
% time and the number of iterations required for convergence for each case
% and visualizes the results.
%
% Author: [Lizhong Jiang]
% Contact: [LizhongJiang_adm@163.com]
% Date: [2025/11/27]
% =========================================================================

%% --- 1. Setup ---
clear;
clc;
close all;

addpath("ProblemFormulation_Constant_VDPO");

fprintf('--- Starting Grid Sensitivity Benchmark for Upper-Bound Method ---\n');

% Define the range of initial interval numbers to test
TMP = 0:10;
initial_interval_counts = 2.^TMP;
initial_interval_counts = [2 5 10 20 40 80 160 320 640 1280];
num_tests = length(initial_interval_counts);

% Pre-allocate arrays to store results
cpu_times = zeros(1, num_tests);
iteration_counts = zeros(1, num_tests);

% Load the base problem and options once
problem_template = define_VDPO_problem();
options = define_solver_options();

% IMPORTANT: Turn off solver display for automated testing to keep the
% command window clean and potentially speed up the process.
options.fmincon.Display = 'none';

%% --- 2. Main Testing Loop ---
for i = 1:num_tests
    
    current_interval_count = initial_interval_counts(i);
    fprintf('\n--- Testing with %d initial interval(s) ---\n', current_interval_count);
    
    % --- Create a copy of the problem and modify it for this test ---
    problem = problem_template;
    
    % Set the initial constraint grid for this specific test case
    problem.pathConstraints.constraintGrid = linspace(problem.time.T0, problem.time.TF, current_interval_count + 1);
    
    % --- Run the solver ---
    try
        [solution, stats] = solve_with_UBMethod(problem, options);
        
        if solution.exitFlag > 0
            % --- Store the results ---
            cpu_times(i) = stats.cpuTime;
            
            % The number of iterations is the number of times the solver was called.
            % This is simply the length of the history arrays.
            iteration_counts(i) = length(stats.history.CnstrGrid);
            
            fprintf('Success! Time: %.3f s, Iterations: %d\n', cpu_times(i), iteration_counts(i));
        else
            % Solver failed to converge
            fprintf('Warning: Solver failed for N = %d (Exit Flag: %d)\n', current_interval_count, solution.exitFlag);
            cpu_times(i) = NaN; % Mark as failed
            iteration_counts(i) = NaN;
        end
        
    catch ME
        % Handle unexpected errors during the solve
        fprintf('ERROR: An error occurred for N = %d: %s\n', current_interval_count, ME.message);
        cpu_times(i) = NaN;
        iteration_counts(i) = NaN;
    end
    
end

fprintf('\n--- Benchmark Complete ---\n');

%% --- 3. Visualization of Results (Combined Dual Y-Axis Plot with Softer Colors) ---
benchmark_grid_visualization

saveResult = false;
if saveResult
    exportgraphics(gcf, "Benchmark_Grid_Ex1_VDPO.pdf", 'Resolution', 600, 'ContentType', 'vector');
end
%% --- 4. Save Workspace Data ---
save('Benchmark_Grid_Results_Ex1_VDPO.mat', 'initial_interval_counts', 'cpu_times', 'iteration_counts');

fprintf('Results saved to Benchmark_Grid_Sensitivity_Results.mat\n');