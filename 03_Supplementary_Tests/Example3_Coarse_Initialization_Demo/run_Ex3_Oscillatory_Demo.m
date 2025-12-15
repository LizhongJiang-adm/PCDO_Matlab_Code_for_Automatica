% =========================================================================
% Script: Oscillatory Constraint Demo (Example 3)
% =========================================================================
% Description:
% This script demonstrates the robustness of the proposed algorithm under
% a coarse initial mesh (Single Interval Initialization).
%
% It solves the PFBF problem (Example 3) starting with the coarsest possible
% constraint partition: Pi^0 = {[T0, TF]}.
%
% This setup typically results in a highly oscillatory path constraint profile
% in the first iteration, verifying the algorithm's ability to handle
% non-monotonic/oscillatory constraints via adaptive refinement.
%
% Author: [Lizhong Jiang]
% Date: [2025/11/17]
% =========================================================================

%% --- 1. Initialization ---
clear; clc; close all;

% Add path to problem-specific functions (Example 3: PFBF)
addpath("ProblemFormulation_Constant_PFBF");

%% --- 2. Setup Problem and Options ---
fprintf('Step 1: Defining the PFBF optimal control problem...\n');
problem = define_PFBF_problem();

fprintf('Step 2: Configuring solver options...\n');
options = define_solver_options();
options.fmincon.Display = 'iter'; % Show iteration details to observe convergence

%% --- 3. Force Coarse Initialization (The Key Step) ---
% Overwrite the default grid to be a single interval [T0, TF]
% This forces the algorithm to start from scratch and handle oscillations.
fprintf('Step 3: Forcing coarse initialization (Single Interval [T0, TF])...\n');
problem.pathConstraints.constraintGrid = [problem.time.T0, problem.time.TF];

%% --- 4. Run Solver ---
fprintf('Step 4: Running the Upper-Bound Method...\n');
[solution, stats] = solve_with_UBMethod(problem, options);

%% --- 5. Analyze Results ---
[g_max, ~, ~, ~, ~] = PathCnstrProfileV3(solution.decisionVariables, problem);
N_CnstrGrid = numel(stats.history.CnstrGrid{end}) - 1;

fprintf('\n---------------------------------------------------\n');
fprintf('Final Results for Oscillatory Demo:\n');
fprintf('Objective Value:       %.3f\n', solution.objectiveValue);
fprintf('Max Constraint Value:  %+.3e (Should be <= 0)\n', g_max);
fprintf('Final Mesh Intervals:  %d\n', N_CnstrGrid);
fprintf('---------------------------------------------------\n');

% Optional: Save the figure of the path constraint at 1st iteration 
% exportgraphics(gcf, "Path_Constraint_Profile_Oscillatory.pdf", ...
%    'Resolution', 600, 'ContentType', 'vector');

%% --- 6. Visualization ---
fprintf('Step 5: Visualizing results...\n');

if solution.exitFlag > 0
    % Plot the results using the shared plotting function
    plot_results(solution, problem, stats);
else
    fprintf('Warning: Solver did not converge.\n');
end

fprintf('\n--- Demo Complete ---\n');