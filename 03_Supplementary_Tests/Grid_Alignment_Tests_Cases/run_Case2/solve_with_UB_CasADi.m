function [solution, stats] = solve_with_UB_CasADi(problem, options)
    % --- 1. 初始化 ---
    startTime = tic;
    iterationCount = 1;
    ocp = problem.ocp;
    grid_points = problem.grid_points; % <-- 从 problem 中获取
    
    % 算法状态：上界约束网格
    N_path_constr = length(problem.path_constraints);
    CnstrGrid = problem.path_constraints_grid;
    % for i = 1:N_path_constr
    %     CnstrGrid{i} = problem.path_constraints.Grid{i}; % <-- 使用 grid_points 初始化
    %     % CnstrGrid{i} = [0,1]; % <-- 使用 grid_points 初始化
    % end
    
    w_opt = []; 
    history = struct();
    

    % --- 2. 主迭代循环 ---
    while true
        fprintf('\n--- Iteration: %d ---\n', iterationCount);
        
        % --- 2a. 构建并转录当前子问题 ---
        grid_points = problem.grid_points;
        strategy = SequentialShootingStrategy_UB(ocp, grid_points);

        for i = 1:N_path_constr
            strategy.add_ineq_path_constraints_UB(problem.path_constraints{i}, CnstrGrid{i});
        end
        strategy.transcribe();
        
        nlp = strategy.get_nlp_problem();
        [obj_fcn, grad_obj_fcn, nonlcon_fcn, grad_nonlcon_fcn, ...
         initial_guess, lbw, ubw] = NmlzNLP(nlp);
        
        if ~isempty(w_opt)
            initial_guess = w_opt; % Warm start
        end

        % --- 2b. 求解 NLP 子问题 ---
        [w_opt, fval, exitflag, ~, lambda] = fmincon(...
            @(w) objective_wrapper(w, obj_fcn, grad_obj_fcn), initial_guess, ...
            [], [], [], [], lbw, ubw, ...
            @(w) nonlcon_wrapper(w, nonlcon_fcn, grad_nonlcon_fcn,options.fmincon.ConstraintTolerance), options.fmincon);
            
        history.CnstrGrid{iterationCount} = CnstrGrid;
        history.iterationPoint{iterationCount} = w_opt;
        history.objVal{iterationCount} = fval;

        % =========================================================================
        % --- 2c. 检查收敛条件 (INTEGRATED LOGIC) ---
        % =========================================================================
        
        if exitflag > 0
            % Only check KKT if fmincon found a feasible solution
            
            % Get all necessary components for the KKT check
            [c_all, ceq_all, grad_c_all, grad_ceq_all] = nonlcon_wrapper(w_opt, nonlcon_fcn, grad_nonlcon_fcn,options.fmincon.ConstraintTolerance);
            [~, grad_f] = objective_wrapper(w_opt, obj_fcn, grad_obj_fcn); 
            
            % Call the dedicated KKT checker
            [isConverged, kktValue] = CheckKKT_UB_CasADi(w_opt, grad_f, c_all, grad_c_all, grad_ceq_all, lbw, ubw, strategy, options);
            
            if isConverged
                break;
                fprintf('\nConvergence criteria met and path constraints satisfied. Terminating.\n');
            end
        end


            
        % =========================================================================
    
        % --- 2d. 更新约束网格 ---
        if exitflag > 0
            % --- Case 1: NLP subproblem was solved successfully ---
            % Refine only the intervals corresponding to active upper-bound constraints.
            
            fprintf('NLP solved. Refining active upper-bound constraint intervals...\n');
            
            % Get the full vector of inequality constraints
            [c_all, ~, ~, ~] = nonlcon_wrapper(w_opt, nonlcon_fcn, grad_nonlcon_fcn,options.fmincon.ConstraintTolerance);
            
            % Determine the number of "other" inequality constraints that come
            % before the upper-bound constraints.
            num_ub_constr = strategy.get_num_ub_constraints();
            num_other_ineq = length(c_all) - num_ub_constr;
            
            offset = num_other_ineq; % Start counting from the end of other constraints
            
            for i = 1:numel(CnstrGrid)
                num_intervals_in_grid = length(CnstrGrid{i}) - 1;
                
                % Extract the upper-bound constraint values for the i-th path constraint
                c_ub_for_path_i = c_all(offset + 1 : offset + num_intervals_in_grid);
                
                % Update the offset for the next loop
                offset = offset + num_intervals_in_grid;
                
                % Find indices of active intervals within this grid
                indices_to_divide = find(c_ub_for_path_i >= -options.fmincon.ConstraintTolerance);
                
                if ~isempty(indices_to_divide)
                    % Calculate midpoints ONLY for the intervals to be divided
                    midpoints_to_add = 0.5 * (CnstrGrid{i}(indices_to_divide) + CnstrGrid{i}(indices_to_divide + 1));
                    
                    % Update the grid
                    CnstrGrid{i} = unique([CnstrGrid{i}, midpoints_to_add]);
                end
            end
            
        else
            % --- Case 2: NLP solver failed (e.g., infeasible) ---
            % Refine all intervals for all path constraints to improve the approximation.
            
            fprintf('NLP failed. Refining all constraint intervals by bisection...\n');
            
            for i = 1:numel(CnstrGrid)
                midpoints = 0.5 * (CnstrGrid{i}(1:end-1) + CnstrGrid{i}(2:end));
                CnstrGrid{i} = unique([CnstrGrid{i}, midpoints]);
            end
        end
        
        iterationCount = iterationCount + 1;
        if iterationCount > 50, break; end
    end
    
    % --- 3. 打包结果 ---
    solution.optimal_weights = w_opt;
    solution.objective_value = fval;

    
    stats.history = history;    
    stats.total_iterations = iterationCount;
    stats.cpu_time = toc(startTime);
end