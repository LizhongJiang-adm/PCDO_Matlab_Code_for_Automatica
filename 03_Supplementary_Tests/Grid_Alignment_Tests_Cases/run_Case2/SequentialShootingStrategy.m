% SequentialShootingStrategy.m (Final version for older CasADi API)
classdef SequentialShootingStrategy < TranscriptionStrategy
    % ... (properties and constructor are unchanged) ...
    properties (SetAccess = protected)
        s_grid
        is_min_time_problem = false;
        Dec
        dae_nmlz
        Ode_reltol = 1e-8;
        Ode_abstol = 1e-8;
    end
    
    methods
        function obj = SequentialShootingStrategy(ocp, grid_points)
            obj@TranscriptionStrategy(ocp, grid_points);
            
            t0 = obj.ocp.t0;
            tf = obj.ocp.tf;
            tol = 1e-9;
            
            if isscalar(tf)
                if abs(grid_points(1) - t0) > tol || abs(grid_points(end) - tf) > tol
                    error(['对于固定时间问题, grid_points 的端点 [%g, %g] 必须与 OCP ' ...
                           '中定义的时间范围 [%g, %g] 精确匹配。'], ...
                           grid_points(1), grid_points(end), t0, tf);
                end
                total_duration = tf - t0;
                obj.s_grid = (grid_points - t0) / total_duration;
                obj.is_min_time_problem = false;
            elseif isvector(tf) && length(tf) == 2
                if abs(grid_points(1) - t0) > tol
                    warning('对于最小时间问题, grid_points(1) 应等于 t0。将以 t0=%g 为准。', t0);
                end
                nominal_duration = grid_points(end) - grid_points(1);
                obj.s_grid = (grid_points - grid_points(1)) / nominal_duration;
                obj.is_min_time_problem = true;
            else
                error('OCP中的tf属性类型无效。它必须是标量或1x2向量。');
            end
        end

        function transcribe(obj)
            import casadi.*
            
            % --- 阶段一: 设置增广系统与创建积分器 (Old API compatible) ---
            fprintf('--- 正在转录 (序贯法): 阶段 1 - 设置积分器 (旧API兼容模式) ---\n');
            mdl = obj.ocp.model;

            % 1.1: 增广系统
            x_sym_integrator = mdl.x_sym;
            ode_integrator = mdl.ode_expr;
            is_augmented = false;
            if ~isempty(obj.ocp.lagrange_term_expr)
                is_augmented = true;
                fprintf('  - 检测到Lagrange项，正在创建增广系统。\n');
                x_L_sym = SX.sym('x_L');
                x_sym_integrator = vertcat(x_sym_integrator, x_L_sym);
                ode_integrator = vertcat(ode_integrator, obj.ocp.lagrange_term_expr);
            end

            % 1.2: 将区间时长dt作为一个参数加入系统
            dt_sym = SX.sym('dt_param');
            ode_scaled = ode_integrator * dt_sym; % dx/ds = f(x,u,p) * dt

            % 1.3: 构建参数列表
            % 参数p现在包括：原始控制u, 原始参数p, 和区间时长dt
            integrator_p_sym = vertcat(mdl.u_sym, mdl.p_sym, dt_sym);
            
            % 1.4: 创建积分器
            dae = struct('x', x_sym_integrator, 'p', integrator_p_sym, 'ode', ode_scaled);
            % 关键：在创建时就固定积分时长为1.0 (s 从 0 到 1)
            integrator_opts = struct('tf', 1.0, 'reltol', obj.Ode_reltol, 'abstol', obj.Ode_abstol);
            integrator_func = integrator('F', 'cvodes', dae, integrator_opts);
            fprintf('  - CasADi积分器 (cvodes) 创建成功。\n');

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
            if is_augmented, X_k{1} = vertcat(X_0, MX(0)); else, X_k{1} = X_0; end
            
            total_duration = T_duration_dec;
            if ~obj.is_min_time_problem, total_duration = obj.ocp.tf - obj.ocp.t0; end

            for k = 1:N
                h_k_norm = obj.s_grid(k+1) - obj.s_grid(k);
                dt_k = total_duration * h_k_norm; % 实际物理时长
                
                % 构建参数向量，现在它必须包含dt_k
                p_k_parts = {};
                if nu > 0, p_k_parts{end+1} = U_k{k}; end
                if np > 0, p_k_parts{end+1} = P_dec; end
                p_k_parts{end+1} = dt_k; % 将时长作为最后一个参数传入
                p_k = vertcat(p_k_parts{:});
                
                % 调用积分器 (不再有 'tf' 参数)
                result = integrator_func('x0', X_k{k}, 'p', p_k);
                X_k{k+1} = result.xf;
            end
            fprintf('  - 前向模拟完成。\n');
            
            % --- 阶段四 & 五: 构建目标/约束并扁平化 (此部分逻辑不变) ---
            % ... (The rest of the code from the previous correct version) ...
            fprintf('--- 正在转录 (序贯法): 阶段 4 & 5 - 构建与扁平化 ---\n');
            nlp_objective = 0;
            g = {}; lbg = []; ubg = [];

            if is_augmented, nlp_objective = nlp_objective + X_k{end}(end); end
            
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
                w0_parts{end+1} = obj.grid_points(end) - obj.grid_points(1);
            end

            obj.nlp_vars_sym = vertcat(w_parts{:});
            obj.nlp_vars_lb = vertcat(lbw_parts{:});
            obj.nlp_vars_ub = vertcat(ubw_parts{:});
            obj.nlp_initial_guess = vertcat(w0_parts{:});
            
            obj.nlp_objective = nlp_objective;
            obj.nlp_constraints = vertcat(g{:});
            obj.nlp_constraints_lb = lbg;
            obj.nlp_constraints_ub = ubg;
            
            obj.Dec.X = X_0;
            obj.Dec.U = U_k;
            obj.Dec.P = P_dec;
            obj.Dec.T = T_duration_dec;
            obj.dae_nmlz = struct('x', x_sym_integrator(1:nx_orig), 'p', integrator_p_sym, 'ode', ode_scaled(1:nx_orig));

            fprintf('转录完成。\n');
        end


        % ... (After the transcribe method) ...
        
        function [X_k_opt, U_k_opt, P_opt, T_duration_opt] = unpack_solution(obj, w_opt)
            % unpack_solution - 从扁平化的解向量中提取结构化的轨迹和参数。
            %
            % (最终方案: 按需重建原始物理系统的积分器)
            % 该方法独立于 transcribe，只关注于重建原始物理系统的轨迹，
            % 忽略任何在优化过程中使用的内部状态增广。
            
            import casadi.*
            
            mdl = obj.ocp.model;
            N = length(obj.s_grid) - 1;
            nx_orig = mdl.nx;
            nu = mdl.nu;
            np = mdl.np;
            
            % --- 阶段一: 解码决策变量 ---
            % 这个过程必须与 transcribe 中的扁平化顺序严格对应
            fprintf('--- 正在解包 ---\n');
            offset = 1;
            X_0_opt = w_opt(offset : offset + nx_orig - 1);
            offset = offset + nx_orig;
            
            U_k_opt = cell(1, N);
            if nu > 0
                for k = 1:N
                    U_k_opt{k} = w_opt(offset : offset + nu - 1);
                    offset = offset + nu;
                end
            end
            
            P_opt = [];
            if np > 0
                P_opt = w_opt(offset : offset + np - 1);
                offset = offset + np;
            end
            
            if obj.is_min_time_problem
                T_duration_opt = w_opt(offset);
            else
                T_duration_opt = obj.ocp.tf - obj.ocp.t0;
            end

            % --- 阶段二: 按需重建原始物理系统的积分器 ---
            
            % 2.1: 只使用模型定义的原始动力学
            x_sym_original = mdl.x_sym;
            ode_original = mdl.ode_expr;
            
            % 2.2: 将dt作为参数对ODE进行时间缩放
            dt_sym = SX.sym('dt_param');
            ode_scaled = ode_original * dt_sym;

            % 2.3: 构建参数列表并创建积分器
            integrator_p_sym = vertcat(mdl.u_sym, mdl.p_sym, dt_sym);
            dae = struct('x', x_sym_original, 'p', integrator_p_sym, 'ode', ode_scaled);
            integrator_opts = struct('tf', 1.0, 'reltol', obj.Ode_reltol, 'abstol', obj.Ode_abstol, 'print_stats', false);
            local_integrator_func = integrator('F_unpack', 'cvodes', dae, integrator_opts);

            % --- 阶段三: 重建物理状态轨迹 ---
            
            X_k_opt = cell(1, N + 1);
            X_k_opt{1} = X_0_opt; % 轨迹从解码出的最优初始状态开始
            
            for k = 1:N
                h_k_norm = obj.s_grid(k+1) - obj.s_grid(k);
                dt_k = T_duration_opt * h_k_norm;
                
                % 构建数值参数向量
                p_k_parts = {};
                if nu > 0, p_k_parts{end+1} = U_k_opt{k}; end
                if np > 0, p_k_parts{end+1} = P_opt; end
                p_k_parts{end+1} = dt_k;
                p_k_numeric = vertcat(p_k_parts{:});
                
                % 调用本地创建的积分器
                result = local_integrator_func('x0', X_k_opt{k}, 'p', p_k_numeric);
                X_k_opt{k+1} = full(result.xf);
            end
        end

        % In SequentialShootingStrategy.m

        % ... (After the unpack_solution method) ...
        
        function [x_dense, u_dense] = get_dense_trajectory(obj, w_opt, t_dense)
            % GET_DENSE_TRAJECTORY - 在任意密集时间点上通过精确积分计算轨迹。
            %
            % 该方法为序贯法提供了一个高保真的轨迹生成工具。它通过在每个
            % 控制区间内，动态地构建并调用一个配置了特定输出网格的CVODES
            % 积分器，来实现对任意密集时间点的精确求解。
            %
            % 核心策略:
            %   1. 解码最优解向量 w_opt，获取最优控制序列和总时长。
            %   2. 遍历每一个控制区间 [t_k, t_{k+1}]。
            %   3. 为该区间内所有需要求解的 t_dense 点创建一个归一化的时间网格。
            %   4. 动态创建一个新的 CasADi 积分器，其 'grid' 选项被设置为
            %      上述的归一化网格。
            %   5. 调用该积分器一次，高效地获得该区间内所有密集点的状态值。
            %
            % 输入:
            %   obj     - (SequentialShootingStrategy) 对象实例。
            %   w_opt   - (vector) NLP求解器返回的、已求解的、扁平化的决策变量向量。
            %   t_dense - (vector) 用户指定的、包含密集时间点的真实时间向量。
            %
            % 输出:
            %   x_dense - (matrix) 状态轨迹矩阵，尺寸为 (nx x length(t_dense))。
            %   u_dense - (matrix) 控制轨迹矩阵，尺寸为 (nu x length(t_dense))。
            
            arguments
                obj (1,1) SequentialShootingStrategy
                w_opt (:,1) {mustBeNumeric}
                t_dense (:,:) {mustBeNumeric, mustBeVector}
            end
            
            % 如果 t_dense 是一个列向量，则将其转置为行向量
            if iscolumn(t_dense)
                t_dense = t_dense';
            end

            import casadi.*
           

            % --- 阶段一: 解码最优决策变量 ---
            % 该过程必须与 transcribe 方法中的扁平化顺序严格对应。
            mdl = obj.ocp.model;
            N = length(obj.s_grid) - 1;
            nx = mdl.nx;
            nu = mdl.nu;
            np = mdl.np;
            
            offset = 1;
            X_0_opt = w_opt(offset : offset + nx - 1);
            offset = offset + nx;
            
            U_k_opt = cell(1, N);
            if nu > 0
                for k = 1:N
                    U_k_opt{k} = w_opt(offset : offset + nu - 1);
                    offset = offset + nu;
                end
            end
            
            P_opt = [];
            if np > 0
                P_opt = w_opt(offset : offset + np - 1);
                offset = offset + np;
            end
            
            if obj.is_min_time_problem
                T_duration_opt = w_opt(offset);
            else
                T_duration_opt = obj.ocp.tf - obj.ocp.t0;
            end

            % --- 阶段二: 准备 DAE 结构体以供后续使用 ---
            % 我们将在循环中基于此 DAE 模板重复创建积分器。
            x_sym_original = mdl.x_sym;
            ode_original = mdl.ode_expr;
            dt_sym = SX.sym('dt_param');
            ode_scaled = ode_original * dt_sym; % 对ODE进行时间缩放
            
            integrator_p_sym = vertcat(mdl.u_sym, mdl.p_sym, dt_sym);
            dae_template = struct('x', x_sym_original, 'p', integrator_p_sym, 'ode', ode_scaled);

            % --- 阶段三: 通过逐段积分重建物理状态轨迹 ---
            % 3.1: 初始化并处理 t0 边界情况
            X_k_opt = cell(1, N + 1);
            X_k_opt{1} = X_0_opt;
            
            x_dense = zeros(nx, length(t_dense));
            u_dense = zeros(nu, length(t_dense));
            
            % 将 t_dense 映射到归一化的总时间 [0, 1]
            t_dense_scaled = (t_dense - obj.ocp.t0) ./ T_duration_opt;

            % cvodes 的 grid 输出不包含 t0 时刻的值，因此需手动处理。
            % 如果 t_dense 中恰好包含初始时刻，我们预先填充x_dense和u_dense在初始时刻的值。
            % 由于CVODE积分器的结果不包含t0点时刻的值，因此在构造grid时，我们只考虑区间(s_start,s_end]上的t_dense点
            % 但是如果t_dense中恰好包含t0值，这样会导致确实一个点的信息
            % 因此在执行前向模拟之前，检查t_dense是否包含t0，如果是，则根据初始信息初始化x_dense和u_dense
            is_t0_in_dense = ismember(obj.ocp.t0, t_dense);
            if any(is_t0_in_dense)
                x_dense(:, is_t0_in_dense) = X_0_opt;
                if nu > 0, u_dense(:, is_t0_in_dense) = U_k_opt{1}; end
            end
            
            % 3.2: 遍历每个控制区间
            for k = 1:N
                % a. 获取当前区间的边界和时长
                s_start = obj.s_grid(k);    % 区间归一化起点
                s_end = obj.s_grid(k+1);      % 区间归一化终点
                h_k_normalized = s_end - s_start;
                dt_k = T_duration_opt * h_k_normalized; % 区间物理时长

                % b. 找到所有落在 (s_start, s_end] 区间内的密集时间点
                % 由于result.xf中不包含起始时刻的值，我们只考虑(t_start, t_end]区间上的ActiveTime点作为grid
                % 这样在取出结果时，我们不需要考虑缺失的起始值
                query_indices = find(t_dense_scaled > s_start & t_dense_scaled <= s_end);
                
                t_in_interval_scaled = t_dense_scaled(query_indices);
                
                % c. 为当前区间动态构建积分器
                % 创建一个包含所有目标点的、相对于区间起点 [0, h_k] 的网格
                grid_current_physical = unique([0, (t_in_interval_scaled - s_start) * T_duration_opt, dt_k]);
                % 将物理时间网格归一化到 [0, 1] 以匹配缩放后的ODE
                grid_current_normalized = grid_current_physical / dt_k;
                
                integrator_opts = struct('grid', grid_current_normalized, 'reltol', obj.Ode_reltol, 'abstol', obj.Ode_abstol, 'print_stats', false);
                local_integrator_func = integrator('F_dense_interval', 'cvodes', dae_template, integrator_opts);
                
                % d. 构建参数并调用积分器
                p_k_parts = {};
                if nu > 0, p_k_parts{end+1} = U_k_opt{k}; end
                if np > 0, p_k_parts{end+1} = P_opt; end
                p_k_parts{end+1} = dt_k;
                p_k_numeric = vertcat(p_k_parts{:});
                
                result = local_integrator_func('x0', X_k_opt{k}, 'p', p_k_numeric);
                
                % e. 填充结果并更新下一个区间的起点
                X_k_opt{k+1} = full(result.xf(:,end));
                
                % 根据 query_indices 是否包含区间的终点来决定提取多少数据
                % 由于t_in_interval中可能包含t_end也可能不包含t_end
                % 如果包含t_end则输出完整result.xf，否则输出result.xf(:,1:end-1)
                num_points_to_fill = length(query_indices);
                x_dense(:, query_indices) = full(result.xf(:, 1:num_points_to_fill));
                if nu > 0
                    u_dense(:, query_indices) = repmat(U_k_opt{k}, 1, num_points_to_fill);
                end
            end
            
            fprintf('  - 密集轨迹生成完毕。\n');
        end        

        function [funExprDec, funExprDec_grad] = getFunAtDiscPoint(obj, Expression, DiscTimeNmlz)
            % GETFUNATDISCPOINT - 获取离散点上表达式关于决策变量的函数
            %
            % 该方法通过一个三阶段过程，高效地计算任意表达式在任意离散
            % 时间点上的值（作为关于所有NLP决策变量的函数）及其梯度。
            %
            % 核心策略 (基于Map):
            %   1. 准备阶段: 收集所有唯一需要计算的时间点，并初始化Map容器。
            %   2. 计算与填充阶段: 通过逐段积分，计算出所有唯一时间点上的
            %      状态，并将 (时间->状态) 和 (时间->控制) 的映射关系存入Map。
            %   3. 分配阶段: 遍历输入表达式，从Map中高效查询所需的状态和
            %      控制，然后构建最终的CasADi函数及其梯度。
            %
            % 输入:
            %   obj           - (object) 类实例。
            %   DiscTimeNmlz  - (cell array) 每个元胞包含一个归一化时间点的向量。
            %   Expression    - (cell array) 与 DiscTimeNmlz 对应的CasADi SX表达式。
            
                        % --- R2019b+ 输入参数验证与规范化 ---
            arguments
                obj (1,1) SequentialShootingStrategy
                Expression (1,:) cell 
                DiscTimeNmlz (1,:) cell
            end

            % --- 自定义验证与修正 ---
            % 1. 验证 Expression 元胞的内容
            isSxExpr = cellfun(@(c) isa(c, 'casadi.SX'), Expression);
            if ~all(isSxExpr)
                error('输入参数 ''Expression'' 的所有元胞必须是 casadi.SX 类型。');
            end
            % 2. 验证 DiscTimeNmlz 元胞的内容并修正形状
            isNumericVectorOrEmpty = cellfun(@(c) (isnumeric(c) && isvector(c)) || isempty(c), DiscTimeNmlz);
            if ~all(isNumericVectorOrEmpty) 
                error('输入参数 ''DiscTimeNmlz'' 的所有元胞必须是数值向量。');
            end
            % 3. 确保 DiscTimeNmlz 中的所有向量都是行向量
            for i = 1:length(DiscTimeNmlz)
                assert(all(0<=DiscTimeNmlz{i}&DiscTimeNmlz{i}<=1), "离散点必须是归一化后在[0,1]区间上的时间点");
                if iscolumn(DiscTimeNmlz{i})
                    DiscTimeNmlz{i} = DiscTimeNmlz{i}';
                end
            end
            % 4. 验证维度匹配
            assert(numel(DiscTimeNmlz) == numel(Expression), "表达式元胞数组与离散点元胞数组的维度必须匹配");



            import casadi.*

            % --- 阶段一: 准备工作 ---
            all_unique_times = unique([DiscTimeNmlz{:}]);
            state_at_time_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
            control_at_time_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
            
            initial_time = obj.s_grid(1);
            if ismember(initial_time, all_unique_times)
                state_at_time_map(initial_time) = obj.Dec.X;
                if obj.ocp.model.nu > 0
                    control_at_time_map(initial_time) = obj.Dec.U{1};
                end
            end
            
            % --- 阶段二: 计算与填充Map ---
            if obj.is_min_time_problem
                total_duration = obj.Dec.T;
            else
                total_duration = obj.ocp.tf - obj.ocp.t0;
            end
            
            X_current = obj.Dec.X;
            N = length(obj.s_grid) - 1;

            for k = 1:N
                s_start = obj.s_grid(k);
                s_end = obj.s_grid(k+1);
                
                times_in_interval = all_unique_times(all_unique_times > s_start & all_unique_times <= s_end);
                
                % 动态创建包含所有目标点的积分器
                % 注意: 即使 times_in_interval 为空，此逻辑依然健壮
                grid_in_interval = unique([s_start, times_in_interval, s_end]);
                grid_normalized = (grid_in_interval - s_start) / (s_end - s_start);
                
                integrator_opts = struct('grid', grid_normalized, 'reltol', obj.Ode_reltol, 'abstol', obj.Ode_abstol, 'print_stats', false);
                local_integrator_func = integrator('F_DiscTime', 'cvodes', obj.dae_nmlz, integrator_opts);

                % 积分并更新下一个起点
                h_k_norm = s_end - s_start;
                dt_k = total_duration * h_k_norm;
                p_k_parts = {};
                if ~isempty(obj.Dec.U) % 这个 if 检查 obj.Dec.U 是否为 {}
                    p_k_parts{end+1} = obj.Dec.U{k}; 
                end
                if ~isempty(obj.Dec.P)
                    p_k_parts{end+1} = obj.Dec.P;
                end
                p_k_parts{end+1} = dt_k;
                p_k = vertcat(p_k_parts{:});

                result = local_integrator_func('x0', X_current, 'p', p_k);
                X_current = result.xf(:, end);

                % 将新计算出的结果填充到Map中
                % 注意: 如果 times_in_interval 为空，此循环不会执行
                for i = 1:length(times_in_interval)
                    time_point = times_in_interval(i);
                    grid_index = find(abs(grid_in_interval - time_point) < 1e-9, 1);
                    
                    state_at_time_map(time_point) = result.xf(:, grid_index - 1);
                    if obj.ocp.model.nu > 0
                        control_at_time_map(time_point) = obj.Dec.U{k};
                    end
                end
            end

            % --- 阶段三: 分配结果 ---
            funExprDec = cell(1, numel(Expression));
            funExprDec_grad = cell(1, numel(Expression));
            
            for i = 1:numel(Expression)
                required_times = DiscTimeNmlz{i};
                
                states_cell = values(state_at_time_map, num2cell(required_times));
                X_for_expr = horzcat(states_cell{:});
                
                U_for_expr = [];
                if obj.ocp.model.nu > 0
                    controls_cell = values(control_at_time_map, num2cell(required_times));
                    U_for_expr = horzcat(controls_cell{:});
                end
                
                expr_func = obj.build_function(Expression{i});
                eval_inputs = obj.get_eval_inputs([], X_for_expr, U_for_expr, obj.Dec.P);
                
                final_expr = [];    final_expr_grad = [];
                if ~isempty(DiscTimeNmlz{i})
                    final_expr = expr_func(eval_inputs{:});
                    final_expr_grad = jacobian(final_expr, obj.nlp_vars_sym);
                end
                
                funExprDec{i} = casadi.Function(['PathFun_', num2str(i)], {obj.nlp_vars_sym}, {final_expr});
                funExprDec_grad{i} = casadi.Function(['GradPathFun_', num2str(i)], {obj.nlp_vars_sym}, {final_expr_grad});
            end
        end
        
    end
    
    methods (Access = protected)
        % ... (helper methods build_function, get_eval_inputs are unchanged) ...
        function casadi_func = build_function(obj, expr)
            import casadi.*
            mdl = obj.ocp.model;
            sym_vars = {mdl.t_sym, mdl.x_sym};
            if mdl.nu > 0, sym_vars{end+1} = mdl.u_sym; end
            if mdl.np > 0, sym_vars{end+1} = mdl.p_sym; end
            casadi_func = Function('func', sym_vars, {expr});
        end

        function inputs = get_eval_inputs(obj, t_val, x_val, u_val, p_val)
            mdl = obj.ocp.model;
            inputs = {t_val, x_val};
            if mdl.nu > 0, inputs{end+1} = u_val; end
            if mdl.np > 0, inputs{end+1} = p_val; end
        end
    end
end