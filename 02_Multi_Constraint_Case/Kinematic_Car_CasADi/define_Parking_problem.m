function problem = define_Parking_problem()
    % --- 1. 定义物理模型 ---
    problem.model = Model_Parking();

        % ocp.add_initial_constraint(x_sym, [-12.8; -3; 0; 0; 0]);
        % problem.N = 10; % 控制区间的数量
        % 单约束


    % --- 2. 定义控制离散化 ---
    problem.N = 10; % 控制区间的数量
    problem.grid_points = linspace(0, 1, problem.N + 1); % 名义时间网格

    % --- 3. 构建 OCP 描述 ---
    ocp = OptimalControlProblem(problem.model);
    x_sym = ocp.model.x_sym;
    
    % --- 4. 设置时间、边界和路径约束表达式 ---
    ocp.set_time_horizon(0, [5, 30]);
    
    ocp.add_initial_constraint(x_sym, [-12.8; -3; 0; 0; 0]);
    % ocp.add_terminal_constraint(x_sym(3), 0, 0);
    ocp.add_terminal_constraint(x_sym(4), 0, 0);
    ocp.add_terminal_constraint(x_sym(1)^2 + x_sym(2)^2, -inf, 0.01);
    
    max_speed = 1.5; max_phi = 0.714; max_accel = 1.5; max_omega = 1; max_theta = 0.25;
    ocp.set_variable_bounds('u', 'indices', 1, 'lb', -max_accel, 'ub', max_accel);
    ocp.set_variable_bounds('u', 'indices', 2, 'lb', -max_omega, 'ub', max_omega);
    
    problem.path_constraints = {
        x_sym(4)-max_speed, 
        -x_sym(4)-max_speed, 
        x_sym(5)-max_phi, 
        -x_sym(5)-max_phi,
        x_sym(3)-max_theta, 
        -x_sym(3)-max_theta
    };
    
    problem.ocp = ocp;
end