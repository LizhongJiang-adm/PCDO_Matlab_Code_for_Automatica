function problem = define_PFBF_2constr_problem()
    % --- 1. 定义物理模型 ---
    problem.model = Model_PFBF();

        % ocp.add_initial_constraint(x_sym, [-12.8; -3; 0; 0; 0]);
        % problem.N = 10; % 控制区间的数量
        % 单约束


    % --- 2. 定义控制离散化 ---
    problem.N = 50; % 控制区间的数量
    problem.grid_points = linspace(0, 150, problem.N + 1); % 名义时间网格

    % --- 3. 构建 OCP 描述 ---
    ocp = OptimalControlProblem(problem.model);
    x_sym = ocp.model.x_sym;
    
    % --- 4. 设置时间、边界和路径约束表达式 ---
    ocp.set_time_horizon(0, 150);
    
    ocp.add_initial_constraint(x_sym, [1; 0.2; 0; 250]);

    ocp.set_objective([], -x_sym(3));

    problem.path_constraints = {x_sym(1)-40,x_sym(2)-0.5};


    for i = 1:numel(problem.path_constraints)
        problem.path_constraints_grid{i} = linspace(0,150,problem.N+1); % <-- 使用 grid_points 初始化
    end
    
    problem.ocp = ocp;
end

