% SequentialShootingStrategy.m (Final version for older CasADi API)
classdef SequentialShootingStrategy_UB < SequentialShootingStrategy
    % ... (properties and constructor are unchanged) ...
    properties (SetAccess = private)
        ineq_path_constraints_UB = {};
        UBIntvl = {};
        is_UBmtd_apx_apply = false;
    end
    
    methods
        function obj = SequentialShootingStrategy_UB(ocp, grid_points)
            % CONSTRUCTOR - SequentialShootingStrategy_UBv3 的构造函数
            
            % --- 调用父类的构造函数 ---
            % 使用 obj@ParentClass(arguments...) 语法。
            % 这会将 ocp 和 grid_points 传递给 SequentialShootingStrategy 的
            % 构造函数，让父类完成所有的验证、s_grid 的计算等工作。
            obj@SequentialShootingStrategy(ocp, grid_points);
        end

        function transcribe(obj)

            import casadi.*
            
            % --- 阶段一: 设置增广系统与创建积分器 (Old API compatible) ---
            fprintf('--- 正在转录 (序贯法): 阶段 1 - 设置积分器 (旧API兼容模式) ---\n');
            mdl = obj.ocp.model;

            % 1.1: 如果存在拉格朗日项，则增广系统
            x_sym_integrator = mdl.x_sym;
            ode_integrator = mdl.ode_expr;
            exist_lagrange = false;
            if ~isempty(obj.ocp.lagrange_term_expr)
                exist_lagrange = true;
                fprintf('  - 检测到Lagrange项，正在创建增广系统。\n');
                x_L_sym = SX.sym('x_L');
                x_sym_integrator = vertcat(x_sym_integrator, x_L_sym);
                ode_integrator = vertcat(ode_integrator, obj.ocp.lagrange_term_expr);
                Idx_Lagrange = numel(x_sym_integrator);
            end

            % 如果需要使用上界方法，增广系统状态
            AllGrid = [];
            if obj.is_UBmtd_apx_apply
                % 合并所有需要计算系统状态的时间点
                AllGrid = unique([obj.UBIntvl{:}]);

                % 记录未增广的状态数
                nx_orig_lgl = numel(x_sym_integrator);
                % 记录有多少个通过上界方法处理的路径约束
                n_UBpath = numel(obj.ineq_path_constraints_UB);
                % 增广系统状态
                c_m = 1e-3;
                TMP = [obj.ineq_path_constraints_UB{:}]';
                g_dot_expr = jtimes(TMP,mdl.x_sym,mdl.ode_expr);
                xdot_aug_UB = 0.5 * (sqrt(g_dot_expr.^2 + c_m.^2) + g_dot_expr);
                x_aug_UB = SX.sym('x_UB',numel(obj.ineq_path_constraints_UB),1);
                x_sym_integrator = vertcat(x_sym_integrator, x_aug_UB);
                ode_integrator = vertcat(ode_integrator, xdot_aug_UB);
                Idx_UBmtd = [numel(x_sym_integrator) + 1 - numel(obj.ineq_path_constraints_UB): ...
                    numel(x_sym_integrator)];
            end



            % 1.2: 将区间时长dt作为一个参数加入系统
            dt_sym = SX.sym('dt_param');
            ode_scaled = ode_integrator * dt_sym; % dx/ds = f(x,u,p) * dt

            % 1.3: 构建参数列表
            % 参数p现在包括：原始控制u, 原始参数p, 和区间时长dt
            integrator_p_sym = vertcat(mdl.u_sym, mdl.p_sym, dt_sym);
            dae = struct('x', x_sym_integrator, 'p', integrator_p_sym, 'ode', ode_scaled);
            
            % --- 阶段二: 定义NLP决策变量 ---
            fprintf('--- 正在转录 (序贯法): 阶段 2 - 定义NLP决策变量 ---\n');
            N = length(obj.s_grid) - 1;
            nx_orig = mdl.nx; nu = mdl.nu; np = mdl.np;

            X_0 = MX.sym('X_0', nx_orig);
            U_k = {}; if nu > 0, for k=1:N, U_k{k} = MX.sym(['U_' num2str(k-1)], nu); end, end
            P_dec = []; if np > 0, P_dec = MX.sym('P', np); end
            T_duration_dec = []; 
            if obj.is_min_time_problem, T_duration_dec = MX.sym('T_duration'); end

            % --- 阶段三: 执行前向模拟循环 ---
            fprintf('--- 正在转录 (序贯法): 阶段 3 - 执行前向模拟 ---\n');
            X_k = cell(1, N + 1);
            if exist_lagrange, X_k{1} = vertcat(X_0, MX(0)); else, X_k{1} = X_0; end
            
            % 定义系统运行时间
            if ~obj.is_min_time_problem
                total_duration = obj.ocp.tf - obj.ocp.t0; 
            else
                total_duration = T_duration_dec;
            end

            % 设定系统初始状态
            if obj.is_UBmtd_apx_apply
                X0 = vertcat(X_k{1},MX(n_UBpath,1));
            else
                X0 = X_k{1};
            end
            
            % 如果约束区间包含初始点，则将系统初始状态传入X_Grid，否则初始化X_Grid为空
            for i = 1:numel(obj.ineq_path_constraints_UB)
                if ismember(obj.s_grid(1),obj.UBIntvl{i}) 
                    X_Grid{i} = X0;
                else
                    X_Grid{i} = MX(numel(ode_integrator),0);
                end
            end

            for k = 1:N
                t0_Intv = obj.s_grid(k);  tf_Intv = obj.s_grid(k+1);
                h_k_norm = tf_Intv - t0_Intv;
                dt_k = total_duration * h_k_norm; % 实际物理时长

                p_k_parts = {};
                if nu > 0, p_k_parts{end+1} = U_k{k}; end
                if np > 0, p_k_parts{end+1} = P_dec; end
                p_k_parts{end+1} = dt_k; % 将时长作为最后一个参数传入
                p_k = vertcat(p_k_parts{:});
                
                grid_curnt = [t0_Intv, AllGrid(t0_Intv<AllGrid & AllGrid<=tf_Intv),  tf_Intv];
                grid_curnt = unique(grid_curnt);
                grid_curnt_scl = (grid_curnt-t0_Intv)/(tf_Intv - t0_Intv);


                for i = 1:numel(grid_curnt_scl)-1
                    t0 = grid_curnt_scl(i); tf = grid_curnt_scl(i+1);
                    integrator_opts = struct('t0', t0, 'tf', tf, 'reltol', obj.Ode_reltol, 'abstol', obj.Ode_abstol);
                    integrator_func = integrator('F', 'cvodes', dae, integrator_opts);

                    result = integrator_func('x0', X0, 'p', p_k);

                    X0 = result.xf(:,end);
                    for j = 1:numel(obj.ineq_path_constraints_UB)
                        if ismember(grid_curnt(i+1),obj.UBIntvl{j})
                            X_Grid{j} = [X_Grid{j},result.xf];
                            Idxeff = Idx_UBmtd(j);
                            X0(Idxeff) = 0;
                        end
                    end
                end

                X_k{k+1} = result.xf(:,end);
            end

            fprintf('  - 前向模拟完成。\n');
            
            % --- 阶段四 & 五: 构建目标/约束并扁平化 (此部分逻辑不变) ---
            % ... (The rest of the code from the previous correct version) ...
            fprintf('--- 正在转录 (序贯法): 阶段 4 & 5 - 构建与扁平化 ---\n');
            nlp_objective = 0;
            g = {}; lbg = []; ubg = [];

            if exist_lagrange, nlp_objective = nlp_objective + X_k{end}(Idx_Lagrange); end
            
            if ~isempty(obj.ocp.mayer_term_expr)
                M_func = obj.build_function(obj.ocp.mayer_term_expr);
                eval_inputs = obj.get_eval_inputs([], X_k{end}(1:nx_orig), [], P_dec);
                nlp_objective = nlp_objective + M_func(eval_inputs{:});
            end

            if obj.is_min_time_problem, nlp_objective = nlp_objective + T_duration_dec; end
            
            for i = 1:length(obj.ocp.initial_constraints)
                c = obj.ocp.initial_constraints{i};
                init_func = obj.build_function(c.expr);
                eval_inputs = obj.get_eval_inputs([], X_0, [], P_dec);
                g{end+1} = init_func(eval_inputs{:});
                lbg = [lbg; c.lb]; ubg = [ubg; c.ub];
            end


            for i = 1:length(obj.ineq_path_constraints_UB)
                expr = obj.ineq_path_constraints_UB{i};
                g_func{i} = obj.build_function(expr);

                IdxCtrl = [];
                for j = 1:numel(obj.UBIntvl{i})-1
                    t0_IntvPoly = obj.UBIntvl{i}(j);  tf_IntvPoly = obj.UBIntvl{i}(j+1);
                    % 找出当前区间对应的控制信号
                    TMP = find(obj.s_grid(1:end-1)<=t0_IntvPoly & tf_IntvPoly<=obj.s_grid(2:end));
                    IdxCtrl = [IdxCtrl, TMP];
                    if isempty(TMP)
                        IdxCtrl = [IdxCtrl, 1]; % 如果约束区间横跨了多个控制区间
                    end
                end

                eval_inputs = obj.get_eval_inputs(MX(0),X_Grid{i}(1:nx_orig,1:end-1) , [U_k{IdxCtrl}], P_dec);
                g_intvl_t0 = g_func{i}(eval_inputs{:});

                Idxeff = Idx_UBmtd(i);
                g_UB{i} = (g_intvl_t0 + X_Grid{i}(Idxeff,2:end))';

            end
           
            for i = 1:length(obj.ocp.path_constraints)
                c = obj.ocp.path_constraints{i};
                path_func = obj.build_function(c.expr);
                for k = 2:(N+1)
                    eval_inputs = obj.get_eval_inputs([], X_k{k}(1:nx_orig), [], P_dec);
                    g{end+1} = path_func(eval_inputs{:});
                    lbg = [lbg; c.lb]; ubg = [ubg; c.ub];
                end
            end
            
            for i = 1:length(obj.ocp.terminal_constraints)
                c = obj.ocp.terminal_constraints{i};
                term_func = obj.build_function(c.expr);
                eval_inputs = obj.get_eval_inputs([], X_k{end}(1:nx_orig), [], P_dec);
                g{end+1} = term_func(eval_inputs{:});
                lbg = [lbg; c.lb]; ubg = [ubg; c.ub];
            end


            % -- 为初始状态 X_0 生成高质量初始猜测 --
            fprintf('  - 正在为初始状态 X_0 生成基于初始约束的猜测。\n');
            x0_guess = zeros(nx_orig, 1); % 1. 从一个默认的全零猜测开始
            
            initial_constraints = obj.ocp.initial_constraints;
            x_sym = mdl.x_sym;
            
            % 2. 遍历所有初始约束
            for i = 1:length(initial_constraints)
                c = initial_constraints{i};
                
                % 检查是否为简单的等式约束
                is_equality = isequal(c.lb, c.ub) && all(isfinite(c.lb));
                
                if is_equality
                    % 情况 A: 约束直接作用于整个状态向量 x_sym
                    expr_str = str(c.expr);
                    
                    % 情况 A: 约束直接作用于整个状态向量 x_sym
                    if isequal(expr_str, str(x_sym))
                        x0_guess = c.lb;
                        break; % 找到了最全面的猜测，可以提前退出循环
                    end
                    
                    % 情况 B: 约束作用于单个状态 x_sym(j)
                    for j = 1:nx_orig
                        if isequal(expr_str, str(x_sym(j)))
                            x0_guess(j) = c.lb; % 只更新该状态的猜测值
                        end
                    end
                end
            end
            
            w_parts = {X_0};
            lbw_parts = {obj.ocp.x_lb}; ubw_parts = {obj.ocp.x_ub};
            w0_parts = {x0_guess}; % 3. 使用新生成的 x0_guess

            if nu > 0
                w_parts{end+1} = vertcat(U_k{:});
                lbw_parts{end+1} = repmat(obj.ocp.u_lb, N, 1);
                ubw_parts{end+1} = repmat(obj.ocp.u_ub, N, 1);

                fprintf('  - 正在为控制变量生成基于边界中点的初始猜测。\n');
                u0_template = zeros(nu, 1);
                u_lb = obj.ocp.u_lb;
                u_ub = obj.ocp.u_ub;
                for i = 1:nu
                    if isfinite(u_lb(i)) && isfinite(u_ub(i))
                        u0_template(i) = (u_lb(i) + u_ub(i)) / 2;
                    end
                end
                w0_parts{end+1} = repmat(u0_template, N, 1);
            end
            if np > 0, w_parts{end+1} = P_dec;
                lbw_parts{end+1} = obj.ocp.p_lb; ubw_parts{end+1} = obj.ocp.p_ub;
                w0_parts{end+1} = zeros(np, 1);
            end
            if obj.is_min_time_problem, w_parts{end+1} = T_duration_dec;
                lbw_parts{end+1} = obj.ocp.tf(1) - obj.ocp.t0;
                ubw_parts{end+1} = obj.ocp.tf(2) - obj.ocp.t0;
                w0_parts{end+1} = (obj.ocp.tf(2) + obj.ocp.tf(1))/2;
            end

            obj.nlp_vars_sym = vertcat(w_parts{:});
            obj.nlp_vars_lb = vertcat(lbw_parts{:});
            obj.nlp_vars_ub = vertcat(ubw_parts{:});
            obj.nlp_initial_guess = vertcat(w0_parts{:});
            
            obj.nlp_objective = nlp_objective;

            if obj.is_UBmtd_apx_apply
                obj.nlp_constraints = vertcat(g{:}, g_UB{:});
                N_Poly = numel(vertcat(g_UB{:}));
                obj.nlp_constraints_lb = [lbg ; -Inf*ones(N_Poly,1)];
                obj.nlp_constraints_ub = [ubg ; 0*ones(N_Poly,1)];
            else
                obj.nlp_constraints = vertcat(g{:});
                obj.nlp_constraints_lb = lbg;
                obj.nlp_constraints_ub = ubg;
            end

            obj.Dec.X = X_0;
            obj.Dec.U = U_k;
            obj.Dec.P = P_dec;
            obj.Dec.T = T_duration_dec;
            obj.dae_nmlz = struct('x', x_sym_integrator(1:nx_orig), 'p', integrator_p_sym, 'ode', ode_scaled(1:nx_orig));

            fprintf('转录完成。\n');
        end


        function [PathCnstrDis, GradPathCnstrDis] = getPathCnstrDis(obj, ActiveTimeNmlz)
            for i = 1:numel(obj.ineq_path_constraints_UB)
                [PathCnstrDis(i), GradPathCnstrDis(i)] = obj.getFunAtDiscPoint(obj.ineq_path_constraints_UB(i), ActiveTimeNmlz(i));
            end
        end

        function obj = add_ineq_path_constraints_UB(obj, expr, CnstrIntvl)
            if nargin > 3, error("上界约束方法只支持逐个添加g<=0形式的规范不等式路径约束"); end
            assert(isscalar(expr)&isa(expr,"casadi.SX"),"expr必须是SX表达式。");
            assert(all(CnstrIntvl(2:end)-CnstrIntvl(1:end-1)>=0)&isvector(CnstrIntvl), ...
                "约束区间必须是单调增长向量" )
            
            obj.ineq_path_constraints_UB{end+1} = expr;
            
            if obj.is_min_time_problem
                % 对于最小化时间问题，我们希望输入的区间就是规范化时间区间
                assert(CnstrIntvl(end)<=1&CnstrIntvl(1)>=0, "对于最小化时间问题，约束区间必须在[0,1]上。");
                obj.UBIntvl{end+1} = CnstrIntvl;
            else
                % 对时间区间做规范化处理
                t0 = obj.ocp.t0;
                tf = obj.ocp.tf;
                total_duration = tf - t0;
                obj.UBIntvl{end+1} = (CnstrIntvl - t0) / total_duration;
            end

            obj.is_UBmtd_apx_apply = true;
        end

        function func_g_UB = buildfunc_ineq_path_UB(obj)
            func_g_UB = cell(1,numel(obj.ineq_path_constraints_UB));
            for i = 1:numel(obj.ineq_path_constraints_UB)
                expr = obj.ineq_path_constraints_UB{i};
                func_g_UB{i} = obj.build_function(expr);
            end
        end

        function [PathCnstrMaxVal, PathCnstrMaxTime_nmlz, path_constrs_trajectory, t_trajectory]...
                = getMaxPConIntvls(obj, w_opt)
            N_smp = 100000;
            [~, ~, P_sol, T_sol] = obj.unpack_solution(w_opt);
            t_trajectory = linspace(0, T_sol, N_smp);
            [x_dense, u_dense] = obj.get_dense_trajectory(w_opt, t_trajectory);

            t_dense_nmlz = linspace(0, 1, N_smp);

            func_g_UB = obj.buildfunc_ineq_path_UB();
            eval_inputs = obj.get_eval_inputs(t_dense_nmlz, x_dense, u_dense, P_sol);
            PathCnstrMaxTime_nmlz = cell(1,numel(obj.ineq_path_constraints_UB));
            PathCnstrMaxVal = cell(1,numel(obj.ineq_path_constraints_UB));
            path_constrs_trajectory = cell(1,numel(obj.ineq_path_constraints_UB));

            for i = 1:numel(obj.ineq_path_constraints_UB)
                % 计算稠密时间点上的路径约束轨迹
                path_constrs_trajectory{i} = full(func_g_UB{i}(eval_inputs{:}));
                for j = 1:numel(obj.UBIntvl{i})-1
                    t0_c = obj.UBIntvl{i}(j);    tf_c = obj.UBIntvl{i}(j+1);
                    t_dense_c = t_dense_nmlz(t0_c<=t_dense_nmlz&t_dense_nmlz<=tf_c);
                    pc_trajectory_c = path_constrs_trajectory{i}(t0_c<=t_dense_nmlz&t_dense_nmlz<=tf_c);
                    [MaxVal,Idx0] = max(pc_trajectory_c);
                    PathCnstrMaxTime_nmlz{i} = [PathCnstrMaxTime_nmlz{i}, t_dense_c(Idx0)];%*T_sol + obj.ocp.t0
                    PathCnstrMaxVal{i} = [PathCnstrMaxVal{i}, MaxVal];
                end
            end
        end

        function num = get_num_ub_constraints(obj)
            num = numel([obj.UBIntvl{:}]) - numel(obj.UBIntvl);
        end
        

    end
    

end