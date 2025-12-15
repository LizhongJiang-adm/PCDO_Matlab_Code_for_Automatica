% =========================================================================
% Main script for solving the Van der Pol Oscillator Optimal Control Problem
% =========================================================================
% Description:
% This script sets up the VDPO problem, configures the solver, runs the
% optimization, and visualizes the results. This is the main entry point.
%
% Author: [Lizhong Jiang]
% Contact: [LizhongJiang_adm@163.com]
% Date: [2025/10/29]
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
clearvars -except ExecutionMode initialCnstrGrid;
clc;
close all;

% Add path to problem-specific functions
addpath("ProblemFormulation_Constant_VDPO");

%% --- 3. Setup Problem and Options ---
fprintf('Step 1: Defining the VDPO optimal control problem...\n');
problem = define_VDPO_problem();


fprintf('Step 2: Configuring solver options...\n');
options = define_solver_options();

%% --- 4. Execute Based on Selected Mode ---
switch ExecutionMode
    case 'Solve'
        % --- SINGLE RUN MODE ---
        fprintf('\n--- Execution Mode: Solve ---\n');
        
        % Ensure solver display is on for detailed feedback
        options.fmincon.Display = 'iter';
        
        fprintf('\nStep 3: Running the solver once...\n');
    
        [solution, stats] = solve_with_UBMethod(problem, options);
        
        [g_max, ~, ~, ~, ~] = PathCnstrProfileV3(solution.decisionVariables, problem);
        N_CnstrGrid = numel(stats.history.CnstrGrid{end})-1;
        fprintf('Final objective function value: %.3f\n', solution.objectiveValue);
        fprintf('Path Constraint: Maximum value = %+.3e\n', g_max);
        fprintf('Final number of inequality constraint : %d\n', N_CnstrGrid);
        
        fprintf('\nStep 4: Visualizing results...\n');
        if solution.exitFlag > 0
            plot_results(solution, problem, stats);
            saveResult = false;
            if saveResult
                exportgraphics(gcf, "Path_Constraint_Profile_VDPO.pdf", 'Resolution', 600, 'ContentType', 'vector');
            end
        else
            fprintf('Solver did not converge. Skipping plot generation.\n');
        end
        
    case 'Benchmark'
        % --- PERFORMANCE BENCHMARK MODE ---
        fprintf('\n--- Execution Mode: Benchmark ---\n');
        
        % Ensure solver display is off for accurate timing
        options.fmincon.Display = 'none';
        
        fprintf('\nStep 3: Measuring performance with timeit...\n');
        
        f_to_time = @() solve_with_UBMethod(problem, options);
        execution_time = timeit(f_to_time);
        
        fprintf('---------------------------------------------------\n');
        fprintf('Benchmark Result for VDPO Problem:\n');
        fprintf('Average execution time: %.4f seconds\n', execution_time);
        fprintf('---------------------------------------------------\n');
        
    otherwise
        error("Invalid ExecutionMode selected. Choose 'Solve' or 'Benchmark'.");
end

fprintf('\n--- Workflow Complete ---\n');

if strcmp(ExecutionMode,'Benchmark') 
    save("Result_VDPO_UB");
end