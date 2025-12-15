% =========================================================================
% Main script for solving the JCB Optimal Control Problem
% =========================================================================
% Author: [Lizhong Jiang]
% Contact: [LizhongJiang_adm@163.com]
% Date: [2025/12/9]
% =========================================================================


clear 
clc;
close all;

addpath("ProblemFormulation_Constant_JCB");


fprintf('Step 1: Defining the JCB problem with 1 control and 1 path constraint...\n');
problem = define_JCB_problem();

fprintf('Step 2: Configuring solver options for the Upper-Bound Method...\n');
options = define_solver_options();


% --- 1. Display Initial Grid Configuration ---
fprintf('--- Initial Grid Setup ---\n');
fprintf('Number of Control Intervals   : %d\n', length(problem.control.grid) - 1);
disp('Control Grid (CVP):');
disp(problem.control.grid);

fprintf('Number of Constraint Intervals: %d\n', length(problem.pathConstraints.constraintGrid) - 1);
disp('Initial Constraint Grid:');
disp(problem.pathConstraints.constraintGrid);

% --- 2. Run the Solver ---
fprintf('\n--- Solving... ---\n');
[solution, stats] = solve_with_UBMethod(problem, options);

% --- 3. Display Results ---
if solution.exitFlag > 0
    iterations = length(stats.history.CnstrGrid);
    
    fprintf('\n--- Results for Case 1 ---\n');
    fprintf('Solver Converged Successfully.\n');
    fprintf('Total Iterations Required: %d\n', iterations);
    fprintf('Final Objective Value:     %.4f\n', solution.objectiveValue);
    
    % Check final constraint violation
    [g_max, ~] = PathCnstrProfileV3(solution.decisionVariables, problem);
    plot_results(solution, problem, stats);
    fprintf('Final Max Constraint Value: %+.4e\n', max(g_max));
else
    fprintf('\n--- Results for Case 1 ---\n');
    fprintf('ERROR: Solver did not converge (Exit Flag: %d).\n', solution.exitFlag);
end


