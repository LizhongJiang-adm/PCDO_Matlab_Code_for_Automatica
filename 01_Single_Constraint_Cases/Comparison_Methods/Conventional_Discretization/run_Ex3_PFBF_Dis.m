% =========================================================================
% Main script for solving the PFBF Problem with the Direct Discretization Method
% =========================================================================
% Description:
% This script sets up the PFBF problem, configures the direct discretization
% solver, runs the optimization, and visualizes the results. This serves as
% a benchmark method.
%
% Author: [Lizhong Jiang]
% Contact: [LizhongJiang adm@163.com]
% Date: [2025/11/05]
% =========================================================================

%% --- 1. Configuration ---
ExecutionMode = 'Solve'; % Options: 'Solve', 'Benchmark'
% ExecutionMode = 'Benchmark';

%% --- 2. Initialization ---
clearvars -except ExecutionMode;
clc;
close all;
addpath("ProblemFormulation_Dis_PFBF");


%% --- 3. Setup Problem and Options ---
fprintf('Step 1: Defining the PFBF problem for the Direct Discretization Method...\n');
problem = define_PFBF_problem_Dis();

fprintf('Step 2: Configuring solver options for the Direct Discretization Method...\n');
options = define_solver_options_Dis();


%% --- 4. Execute Based on Selected Mode ---
switch ExecutionMode
    case 'Solve'
        % --- SINGLE RUN MODE ---
        fprintf('\n--- Execution Mode: Solve ---\n');
        options.fmincon.Display = 'iter';
        
        fprintf('\nStep 3: Running the Direct Discretization solver once...\n');
        [solution, stats] = solve_with_DisMethod(problem, options);

        [g_max, ~, ~, ~, ~] = PathCnstrProfileV3(solution.decisionVariables, problem);
        N_CnstrGrid = numel(problem.control.grid);
        fprintf('Final objective function value: %.3f\n', solution.objectiveValue);
        fprintf('Path Constraint: Maximum value = %+.3e\n', g_max);
        fprintf('Final number of inequality constraint : %d\n', N_CnstrGrid);
        
        fprintf('\nStep 4: Visualizing results...\n');
        if solution.exitFlag > 0
            plot_results(solution, problem, stats);
        else
            fprintf('Solver did not converge. Skipping plot generation.\n');
        end
        
    case 'Benchmark'
        % --- PERFORMANCE BENCHMARK MODE ---
        fprintf('\n--- Execution Mode: Benchmark ---\n');
        options.fmincon.Display = 'none';
        
        fprintf('\nStep 3: Measuring performance with timeit...\n');
        f_to_time = @() solve_with_DisMethod(problem, options);
        execution_time = timeit(f_to_time);
        
        fprintf('---------------------------------------------------\n');
        fprintf('Benchmark Result for PFBF (Discretization Method):\n');
        fprintf('Control Intervals (N): %d\n', problem.control.N);
        fprintf('Average execution time: %.4f seconds\n', execution_time);
        fprintf('---------------------------------------------------\n');
        
    otherwise
        error("Invalid ExecutionMode selected. Choose 'Solve' or 'Benchmark'.");
end

fprintf('\n--- PFBF (Discretization Method) Workflow Complete ---\n');

if strcmp(ExecutionMode,'Benchmark') 
    save("Result_PFBF_Dis");
end