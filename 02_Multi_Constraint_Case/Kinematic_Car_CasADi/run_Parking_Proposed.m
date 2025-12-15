% =========================================================================
% Main script for solving the Parking Problem with the Upper-Bound Method
% =========================================================================
% Description:
% This script sets up the Parking Problem, configures the Upper-Bound
% solver, runs the optimization, and visualizes the results.
% It supports 'Solve' and 'Benchmark' execution modes.
%
% Author: [Lizhong Jiang]
% Contact: [LizhongJiang_adm@163.com]
% Date: [2025/11/13]
% =========================================================================

%% --- 1. Configuration ---
ExecutionMode = 'Solve'; % Options: 'Solve', 'Benchmark'
% ExecutionMode = 'Benchmark'; % Options: 'Solve', 'Benchmark'

%% --- 2. Initialization ---
clearvars -except ExecutionMode;
clc;
close all;



try
    import casadi.*
catch
    error('CasADi not found! Please download and add CasADi to your MATLAB path.');
end


%% --- 3. Setup Problem and Options ---
fprintf('Step 1: Defining the Parking problem...\n');
problem = define_Parking_problem();

fprintf('Step 2: Configuring solver options for the Upper-Bound Method...\n');
options = define_solver_options_UB();


%% --- 4. Execute Based on Selected Mode ---
switch ExecutionMode
    case 'Solve'
        % --- SINGLE RUN MODE with VISUALIZATION ---
        fprintf('\n--- Execution Mode: Solve ---\n');
        options.fmincon.Display = 'iter';
        
        fprintf('Step 3: Running the Upper-Bound solver...\n');
        [solution, stats] = solve_with_UB_CasADi(problem, options);
        
        fprintf('\nStep 4: Visualizing results...\n');

        % --- Visualization Logic ---
        % Recreate the strategy object to access helper/visualization methods.
        % This is necessary because the final 'CnstrGrid' is inside the solver's stats.
        
        w_opt = solution.optimal_weights;
        ocp = problem.ocp;
        grid_points = problem.grid_points;
        
        % Use the final constraint grid from the solver's history for full accuracy
        final_CnstrGrid = stats.history.CnstrGrid{end};
        
        strategy = SequentialShootingStrategy_UB(ocp, grid_points);
        % Re-add constraints with the FINAL grid to have a complete strategy object
        for i = 1:numel(problem.path_constraints)
            strategy.add_ineq_path_constraints_UB(problem.path_constraints{i}, final_CnstrGrid{i});
        end
        strategy.transcribe();
        
        saveResult = false; 
        visualization_ParkingV2
        visualization_Iteration

    
    case 'Benchmark'
        % --- PERFORMANCE BENCHMARK MODE ---
        fprintf('\n--- Execution Mode: Benchmark ---\n');
        options.fmincon.Display = 'none';
        
        fprintf('Step 3: Measuring performance with timeit...\n');
        f_to_time = @() solve_with_UB_CasADi(problem, options);
        execution_time = timeit(f_to_time);
        
        fprintf('---------------------------------------------------\n');
        fprintf('Benchmark Result for Parking (Upper-Bound Method):\n');
        fprintf('Control Intervals (N): %d\n', problem.N);
        fprintf('Average execution time: %.4f seconds\n', execution_time);
        fprintf('---------------------------------------------------\n');
        
    otherwise
        error("Invalid ExecutionMode selected. Choose 'Solve' or 'Benchmark'.");
end

fprintf('\n--- Parking (Upper-Bound Method) Workflow Complete ---\n');

if strcmp(ExecutionMode,'Benchmark') 
    save("Result_Parking_UB");
end