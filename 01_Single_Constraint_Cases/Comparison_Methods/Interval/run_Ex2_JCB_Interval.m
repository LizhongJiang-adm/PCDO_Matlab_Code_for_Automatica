% =========================================================================
% Main script for solving the JCB Problem with the Interval Method
% =========================================================================
% Description:
% This script sets up the JCB problem, configures the Interval Method solver,
% runs the optimization, and visualizes the results. It supports multiple
% execution modes for either solving the problem or benchmarking performance.
%
% Author: [Lizhong Jiang]
% Contact: [LizhongJiang adm@163.com]
% Date: [2025/10/30]
% =========================================================================

%% --- 1. Configuration ---
% =========================================================================
% Set the execution mode for this script.
% 'Solve'     - Runs the solver once to get a solution and generate plots.
% 'Benchmark' - Runs the solver multiple times using 'timeit' for an
%               accurate performance measurement. No plots are generated.
% =========================================================================
ExecutionMode = 'Solve'; % Options: 'Solve', 'Benchmark'
% ExecutionMode = 'Benchmark'; 

%% --- 2. Initialization ---
clearvars -except ExecutionMode;
clc;
close all;

% Add path to the folder containing the JCB problem-specific functions
addpath("ProblemFormulation_Intvl_JCB");


%% --- 3. Setup Problem and Options ---
fprintf('Step 1: Defining the JCB problem for the Interval Method...\n');
problem = define_JCB_problem_Interval();

fprintf('Step 2: Configuring solver options for the Interval Method...\n');
options = define_solver_options_Interval();


%% --- 4. Execute Based on Selected Mode ---
switch ExecutionMode
    case 'Solve'
        % --- SINGLE RUN MODE ---
        fprintf('\n--- Execution Mode: Solve ---\n');
        
        % Ensure solver display is on for detailed feedback
        options.fmincon.Display = 'iter';
        
        fprintf('\nStep 3: Running the Interval Method solver once...\n');
        [solution, stats] = solve_with_IntervalMethod(problem, options);

        [g_max, ~, ~, ~, ~] = PathCnstrProfileV3(solution.decisionVariables, problem);
        N_CnstrGrid = numel(stats.history.constraintGrid{end})-1;
        fprintf('Final objective function value: %.3f\n', solution.objectiveValue);
        fprintf('Path Constraint: Maximum value = %+.3e\n', g_max);
        fprintf('Final number of inequality constraint : %d\n', N_CnstrGrid);

        fprintf('\nStep 4: Visualizing results...\n');
        if solution.exitFlag > 0
            % Reuse the standardized plotting function
            plot_results(solution, problem, stats);
        else
            fprintf('Solver did not converge. Skipping plot generation.\n');
        end
        
    case 'Benchmark'
        % --- PERFORMANCE BENCHMARK MODE ---
        fprintf('\n--- Execution Mode: Benchmark ---\n');
        
        % Ensure solver display is off for accurate timing
        options.fmincon.Display = 'none';
        
        fprintf('\nStep 3: Measuring performance with timeit...\n');
        
        f_to_time = @() solve_with_IntervalMethod(problem, options);
        execution_time = timeit(f_to_time);
        
        fprintf('---------------------------------------------------\n');
        fprintf('Benchmark Result for JCB (Interval Method):\n');
        fprintf('Control Intervals (N): %d\n', problem.control.N);
        fprintf('Average execution time: %.4f seconds\n', execution_time);
        fprintf('---------------------------------------------------\n');
        
    otherwise
        error("Invalid ExecutionMode selected. Choose 'Solve' or 'Benchmark'.");
end

fprintf('\n--- JCB (Interval Method) Workflow Complete ---\n');


if strcmp(ExecutionMode,'Benchmark') 
    save("Result_JCB_Interval");
end