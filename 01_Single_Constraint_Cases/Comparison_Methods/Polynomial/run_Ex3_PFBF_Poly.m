% =========================================================================
% Main script for solving the PFBF Problem with the Polynomial Method
% =========================================================================
% Description:
% This script sets up the PFBF problem, configures the Polynomial
% Approximation solver, runs the optimization, and visualizes the results.
%
% Author: [Lizhong Jiang]
% Contact: [LizhongJiang_adm@163.com]
% Date: [2025/10/30]
% =========================================================================

%% --- 1. Configuration ---
ExecutionMode = 'Solve'; % Options: 'Solve', 'Benchmark'
% ExecutionMode = 'Benchmark';

%% --- 2. Initialization ---
clearvars -except ExecutionMode;
clc;
close all;
addpath("ProblemFormulation_Poly_PFBF") 


%% --- 3. Setup Problem and Options ---
fprintf('Step 1: Defining the PFBF problem for the Polynomial Method...\n');
problem = define_PFBF_problem_Poly();

fprintf('Step 2: Configuring solver options for the Polynomial Method...\n');
options = define_solver_options_Poly();


%% --- 4. Execute Based on Selected Mode ---
switch ExecutionMode
    case 'Solve'
        % --- SINGLE RUN MODE ---
        fprintf('\n--- Execution Mode: Solve ---\n');
        options.fmincon.Display = 'iter';
        
        fprintf('\nStep 3: Running the Polynomial Method solver once...\n');
        [solution, stats] = solve_with_PolyMethod(problem, options);

        [g_max, ~, ~, ~, ~] = PathCnstrProfileV3(solution.decisionVariables, problem);
        N_CnstrGrid = numel(stats.history.T_Dis{end})-1;
        fprintf('Final objective function value: %.3f\n', solution.objectiveValue);
        fprintf('Path Constraint: Maximum value = %+.3e\n', g_max);
        fprintf('Final number of inequality constraint : %d\n', 2*N_CnstrGrid);
        
        fprintf('\nStep 4: Visualizing results...\n');
        if solution.exitFlag > 0
            % A dedicated plotting function for the Polynomial method is recommended
            % to visualize the polynomial approximations as well.
            plot_results(solution, problem, stats);
        else
            fprintf('Solver did not converge. Skipping plot generation.\n');
        end
        
    case 'Benchmark'
        % --- PERFORMANCE BENCHMARK MODE ---
        fprintf('\n--- Execution Mode: Benchmark ---\n');
        options.fmincon.Display = 'none';
        
        fprintf('\nStep 3: Measuring performance with timeit...\n');
        f_to_time = @() solve_with_PolyMethod(problem, options);
        execution_time = timeit(f_to_time);
        
        fprintf('---------------------------------------------------\n');
        fprintf('Benchmark Result for PFBF (Polynomial Method):\n');
        fprintf('Control Intervals (N): %d\n', problem.control.N);
        fprintf('Average execution time: %.4f seconds\n', execution_time);
        fprintf('---------------------------------------------------\n');
        
    otherwise
        error("Invalid ExecutionMode selected. Choose 'Solve' or 'Benchmark'.");
end

fprintf('\n--- PFBF (Polynomial Method) Workflow Complete ---\n');

if strcmp(ExecutionMode,'Benchmark') 
    save("Result_PFBF_Poly");
end