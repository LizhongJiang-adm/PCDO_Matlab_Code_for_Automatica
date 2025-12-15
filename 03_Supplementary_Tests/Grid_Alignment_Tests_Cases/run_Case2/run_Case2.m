% =========================================================================
% Script: PFBF with 2 Path Constraints 
% =========================================================================
%
% Author: [Lizhong Jiang]
% Contact: [LizhongJiang_adm@163.com]
% Date: [2025/11/29]
% =========================================================================

%% --- 1. Initialization ---
clear;
clc;
close all;

try
    import casadi.*
catch
    error('CasADi not found! Please download and add CasADi to your MATLAB path.');
end


%% --- 2. Setup Problem and Options ---
fprintf('Step 1: Defining the PFBF problem with 2 path constraints...\n');
problem = define_PFBF_2constr_problem();

fprintf('Step 2: Configuring solver options for the Upper-Bound Method...\n');
options = define_solver_options_UB();

% Suppress detailed fmincon output to keep the summary clean
options.fmincon.Display = 'iter';


%% --- 3. Display Initial Grid Configuration ---
fprintf('\n======================================================\n');
fprintf('        Initial Grid Configuration for PFBF (2 Constr)\n');
fprintf('======================================================\n');
fprintf('--- Control Grid ---\n');
fprintf('Number of Control Intervals: %d\n', length(problem.grid_points) - 1);
disp('Control Grid Points:');
disp(problem.grid_points);

fprintf('\n--- Constraint Grids ---\n');
fprintf('Number of Constraint 1 Intervals: %d\n', length(problem.path_constraints_grid{1}) - 1);
disp('Initial Constraint Grid 1 Points:');
disp(problem.path_constraints_grid{1});

fprintf('Number of Constraint 2 Intervals: %d\n', length(problem.path_constraints_grid{2}) - 1);
disp('Initial Constraint Grid 2 Points:');
disp(problem.path_constraints_grid{2});

%% --- 4. Run the Solver ---
fprintf('\n--- Solving... ---\n');
[solution, stats] = solve_with_UB_CasADi(problem, options);

%% --- 5. Display Final Results ---

iterations = length(stats.history.CnstrGrid);

fprintf('\n======================================================\n');
fprintf('                      Final Results\n');
fprintf('======================================================\n');
fprintf('Solver Converged Successfully.\n');
fprintf('Total Iterations Required: %d\n', iterations);
fprintf('Final Objective Value:     %.4f\n', solution.objective_value);
