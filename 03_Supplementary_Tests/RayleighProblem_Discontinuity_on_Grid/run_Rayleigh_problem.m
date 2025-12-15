% =========================================================================
% Script: Rayleigh Problem 
% =========================================================================
% Description:
% This script runs the Rayleigh optimal control problem 
%
% Author: [Lizhong Jiang]
% Contact: [LizhongJiang_adm@163.com]
% Date: [Date]
% =========================================================================

%% --- 1. Initialization ---
clear;
clc;
close all;

% Add path to the problem-specific functions
addpath("ProblemFormulation_Constant_RAYL");

%% --- 2. Setup Problem and Options ---
fprintf('Step 1: Defining the Rayleigh optimal control problem...\n');
problem = define_RAYL_problem();

fprintf('Step 2: Configuring solver options for the Upper-Bound Method...\n');
options = define_solver_options();

% Ensure fmincon provides iterative feedback to observe the process
options.fmincon.Display = 'iter';

%% --- 3. Run the Solver ---
fprintf('\nStep 3: Running the Upper-Bound Method on the Rayleigh problem...\n');
[solution, stats] = solve_with_UBMethod(problem, options);

%% --- 4. Display Final Summary ---
if solution.exitFlag > 0
    % Calculate final constraint violation for reporting
    [g_profile_max, t_at_max, state_and_g, time_vector, control_vector] = ...
        PathCnstrProfileV3(solution.decisionVariables, problem);

    final_grid_points = stats.history.CnstrGrid{end};
    
    fprintf('\n-------------------- Demo Summary --------------------\n');
    fprintf('Solver converged successfully.\n');
    fprintf('Final Objective Value:      %.4f\n', solution.objectiveValue);
    fprintf('Maximum Constraint Value:   %+.4e\n', max(g_profile_max));
    fprintf('Total Iterations:           %d\n', stats.totalIterations);
    fprintf('Final Constraint Intervals: %d\n', length(final_grid_points) - 1);
    fprintf('------------------------------------------------------\n');
else
    fprintf('\n--- Solver did not converge. ---\n');
end

%% --- 5. Visualization ---
fprintf('\nStep 5: Visualizing results...\n');

plot_constraint_composition(solution, problem);
%exportgraphics(gcf, "Path_Constraint_Rayleigh_problem.pdf", 'Resolution', 600, 'ContentType', 'vector');

