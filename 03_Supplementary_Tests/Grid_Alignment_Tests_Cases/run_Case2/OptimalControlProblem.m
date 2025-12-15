% OptimalControlProblem.m
classdef OptimalControlProblem < handle
    % OptimalControlProblem - 以标准形式描述一个连续时间的最优控制问题。
    %
    % 这是一个纯粹的问题定义类，它将动力学模型与目标、约束和变量边界相结合。
    
    properties (SetAccess = protected)
        model
        t0 = 0
        % tf - 终端时间。可以是一个标量(固定时间)或一个1x2向量[min, max](最小时间)。
        tf = [] 
        lagrange_term_expr
        mayer_term_expr
        path_constraints = {}
        initial_constraints = {}
        terminal_constraints = {}
        x_lb = [], x_ub = []
        u_lb = [], u_ub = []
        p_lb = [], p_ub = []
    end
    
    methods
        function obj = OptimalControlProblem(model)
            if ~isa(model, 'DynamicModel')
                error('输入模型必须是 DynamicModel 的子类。');
            end
            obj.model = model;
            obj.x_lb = -inf(model.nx, 1); obj.x_ub = inf(model.nx, 1);
            obj.u_lb = -inf(model.nu, 1); obj.u_ub = inf(model.nu, 1);
            obj.p_lb = -inf(model.np, 1); obj.p_ub = inf(model.np, 1);
        end
        
        function obj = set_objective(obj, lagrange_expr, mayer_expr)
            if nargin < 2, lagrange_expr = []; end
            if nargin < 3, mayer_expr = []; end
            obj.lagrange_term_expr = lagrange_expr;
            obj.mayer_term_expr = mayer_expr;
        end

        function obj = add_path_constraint(obj, expr, lb, ub)
            if nargin < 3, lb = -inf; end
            if nargin < 4, ub = 0; end
            obj.path_constraints = obj.add_constraint_internal(obj.path_constraints, expr, lb, ub);
        end
        
        function obj = add_initial_constraint(obj, expr, lb, ub)
            if nargin < 3, lb = 0; end
            if nargin < 4, ub = lb; end
            obj.initial_constraints = obj.add_constraint_internal(obj.initial_constraints, expr, lb, ub);
        end

        function obj = add_terminal_constraint(obj, expr, lb, ub)
            if nargin < 3, lb = 0; end
            if nargin < 4, ub = lb; end
            obj.terminal_constraints = obj.add_constraint_internal(obj.terminal_constraints, expr, lb, ub);
        end

        function obj = set_variable_bounds(obj, type, varargin)
            p = inputParser;
            addRequired(p, 'type', @(x) ismember(x, {'x', 'u', 'p'}));
            addParameter(p, 'indices', 'all', @(x) ischar(x) || isvector(x));
            addParameter(p, 'lb', -inf);
            addParameter(p, 'ub', inf);
            parse(p, type, varargin{:});
            res = p.Results;
            prop_lb = [res.type, '_lb'];
            prop_ub = [res.type, '_ub'];
            if ischar(res.indices) && strcmp(res.indices, 'all')
                indices = 1:obj.model.(['n' res.type]);
            else, indices = res.indices; 
            end
            if ~isinf(res.lb), obj.(prop_lb)(indices) = res.lb; end
            if ~isinf(res.ub), obj.(prop_ub)(indices) = res.ub; end
        end

        function obj = set_time_horizon(obj, t0, tf)
            % set_time_horizon - 设置初始时间和终端时间。
            %
            % 用法:
            %   ocp.set_time_horizon(0, 10)       % 固定时间问题, tf = 10
            %   ocp.set_time_horizon(0, [10, 20]) % 有界最小时间问题, tf in [10, 20]
            %   ocp.set_time_horizon(0)           % 无界最小时间问题, tf >= t0
            
            validateattributes(t0, {'numeric'}, {'scalar', 'real'}, 'set_time_horizon', 't0');
            obj.t0 = t0;
            
            if nargin < 3 || isempty(tf)
                % 无界最小时间问题
                obj.tf = [obj.t0, inf];
            elseif isscalar(tf)
                % 固定终端时间问题
                if tf <= obj.t0
                    error('固定终端时间 tf 必须大于初始时间 t0。');
                end
                obj.tf = tf;
            elseif isvector(tf) && length(tf) == 2
                % 有界最小时间问题
                tf_min = tf(1);
                tf_max = tf(2);
                if tf_min < obj.t0 || tf_max < tf_min
                    error('无效的时间边界 [tf_min, tf_max]。必须满足 t0 <= tf_min <= tf_max。');
                end
                obj.tf = [tf_min, tf_max];
            else
                error('无效的 tf 参数。它必须是标量、[min, max] 向量或空值。');
            end
        end

        function [eq_exprs, ineq_exprs] = normalize_constraints(~, constraints_cell)
            eq_exprs = {}; ineq_exprs = {};
            for i = 1:numel(constraints_cell)
                g = constraints_cell{i}.expr; lb = constraints_cell{i}.lb; ub = constraints_cell{i}.ub;
                for j = 1:numel(g)
                    g_j = g(j); lb_j = lb(j); ub_j = ub(j);
                    if lb_j == ub_j
                        eq_exprs{end+1} = g_j - lb_j;
                    else
                        if ~isinf(ub_j), ineq_exprs{end+1} = g_j - ub_j; end
                        if ~isinf(lb_j), ineq_exprs{end+1} = lb_j - g_j; end
                    end
                end
            end
        end

    end

    methods (Access = private)
        function cell_array = add_constraint_internal(~, cell_array, expr, lb, ub)
            expr_size = size(expr);
            if isscalar(lb), lb = lb * ones(expr_size);
            elseif ~isequal(size(lb), expr_size), error('约束表达式的维度与下界的维度不匹配。'); end
            if isscalar(ub), ub = ub * ones(expr_size);
            elseif ~isequal(size(ub), expr_size), error('约束表达式的维度与上界的维度不匹配。'); end
            constraint.expr = expr;
            constraint.lb = lb;
            constraint.ub = ub;
            cell_array{end+1} = constraint;
        end

        
    end
end