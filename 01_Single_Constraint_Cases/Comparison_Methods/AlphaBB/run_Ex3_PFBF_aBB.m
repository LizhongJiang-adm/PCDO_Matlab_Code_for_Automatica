% =========================================================================
% Main script for solving the PFBF Problem with the aBB Method
% =========================================================================
% Description:
% This script sets up the PFBF problem, configures the aBB Method solver,
% runs the optimization, and visualizes the results.
%
% Author: [Lizhong Jiang]
% Contact: [LizhongJiang_adm@163.com]
% Date: [2025/11/8]
% =========================================================================

%% --- 1. Configuration ---
ExecutionMode = 'Solve'; % Options: 'Solve', 'Benchmark'
% ExecutionMode = 'Benchmark';

%% --- 2. Initialization ---
clearvars -except ExecutionMode;
clc;
close all;
addpath("ProblemFormulation_aBB_PFBF");


%% --- 3. Setup Problem and Options ---
fprintf('Step 1: Defining the PFBF problem for the aBB Method...\n');
problem = define_PFBF_problem_aBB();

fprintf('Step 2: Configuring solver options for the aBB Method...\n');
options = define_solver_options_aBB();


%% --- 4. Execute Based on Selected Mode ---
switch ExecutionMode
    case 'Solve'
        % --- SINGLE RUN MODE ---
        fprintf('\n--- Execution Mode: Solve ---\n');
        options.fmincon.Display = 'iter';
        
        fprintf('\nStep 3: Running the aBB Method solver once...\n');
        [solution, stats] = solve_with_aBBMethod(problem, options);
        
        [g_max, ~, ~, ~, ~] = PathCnstrProfileV3(solution.decisionVariables, problem);
        N_CnstrGrid = numel(stats.history.eta{end})-1;
        fprintf('Final objective function value: %.3f\n', solution.objectiveValue);
        fprintf('Path Constraint: Maximum value = %+.3e\n', g_max);
        fprintf('Final number of inequality constraint : %d\n', N_CnstrGrid);
        fprintf('Final number of equality constraint : %d\n', 3*N_CnstrGrid);
        fprintf('Final number of decision variables : %d\n', problem.control.N + 3*N_CnstrGrid);

        fprintf('\nStep 4: Visualizing results...\n');
        if solution.exitFlag > 0
            % Reuse the standardized plotting function.
            plot_results(solution, problem, stats);
        else
            fprintf('Solver did not converge. Skipping plot generation.\n');
        end
        
    case 'Benchmark'
        % --- PERFORMANCE BENCHMARK MODE ---
        fprintf('\n--- Execution Mode: Benchmark ---\n');
        options.fmincon.Display = 'none';
        
        fprintf('\nStep 3: Measuring performance with timeit...\n');
        f_to_time = @() solve_with_aBBMethod(problem, options);
        execution_time = timeit(f_to_time);
        
        fprintf('---------------------------------------------------\n');
        fprintf('Benchmark Result for PFBF (aBB Method):\n');
        fprintf('Control Intervals (N): %d\n', problem.control.N);
        fprintf('Average execution time: %.4f seconds\n', execution_time);
        fprintf('---------------------------------------------------\n');
        
    otherwise
        error("Invalid ExecutionMode selected. Choose 'Solve' or 'Benchmark'.");
end

fprintf('\n--- PFBF (aBB Method) Workflow Complete ---\n');

if strcmp(ExecutionMode,'Benchmark') 
    save("Result_PFBF_aBB");
end