% DynamicModel.m
classdef (Abstract) DynamicModel < handle
    % DynamicModel - 定义一个通用、健壮的动力学模型的抽象基类
    %
    % 这是一个“契约”，所有具体的物理模型（如ParkingModel, VdpModel）
    % 都应该继承此类，并实现其抽象方法。它被设计为一个句柄类，
    % 以便在程序中高效地传递引用，并允许方法修改对象状态。
    %
    % 核心功能:
    % 1. 定义了一套标准的属性来存储模型的符号定义。
    % 2. 强制子类通过实现 setup_model 方法来提供具体的模型方程。
    % 3. 自动计算模型维度 (nx, nu, np)。
    % 4. 提供一个便利方法 get_dynamics_function 来获取一个可调用的
    %    、高效的 CasADi 函数，该函数签名会智能地适应模型是否
    %    包含控制(u)或参数(p)。
    
    properties (SetAccess = protected)
        % --- 符号变量定义 ---
        x_sym       % 状态向量 (CasADi SX)
        u_sym       % 控制向量 (CasADi SX) - 可以为空
        p_sym       % 模型参数 (CasADi SX, 向量或结构体) - 可以为空
        t_sym       % 时间变量 (CasADi SX)
        
        % --- 符号表达式 ---
        ode_expr    % 核心ODE表达式: dx/dt = f(t, x, u, p)
        
        % --- 维度信息 (由构造函数自动计算) ---
        nx = 0      % 状态数量
        nu = 0      % 控制数量
        np = 0      % 参数数量
    end
    
    methods (Abstract)
        % --- 抽象方法：子类的强制性任务 ---
        
        % setup_model - 定义模型的具体内容
        %
        % 这是所有子类都必须实现的“契约”方法。在这个方法中，
        % 子类必须为 x_sym, u_sym, p_sym, t_sym, 和 ode_expr 
        % 这些核心属性赋上具体的值。
        setup_model(obj) 
    end
    
    methods
        % --- 具体方法：父类提供的通用功能 ---
        
        function obj = DynamicModel()
            % 构造函数 (Constructor)
            import casadi.*
            
            % 1. 为核心属性设置安全的默认值
            obj.x_sym = [];
            obj.u_sym = [];
            obj.p_sym = [];
            obj.t_sym = SX.sym('t'); % 提供一个默认的时间符号变量
            obj.ode_expr = [];
            
            % 2. 调用子类实现的 setup_model 方法来填充模型定义
            obj.setup_model();
            
            % 3. 在模型定义完成后，健壮地计算并存储维度信息
            if ~isempty(obj.x_sym)
                obj.nx = size(obj.x_sym, 1);
            end
            if ~isempty(obj.u_sym)
                obj.nu = size(obj.u_sym, 1);
            end
            if ~isempty(obj.p_sym)
                obj.np = size(obj.p_sym, 1);
            end
        end
        
        function dyn_fun = get_dynamics_function(obj, func_name)
            % GET_DYNAMICS_FUNCTION - 将符号表达式编译成一个高效的CasADi函数
            %
            % 这个函数返回一个可直接调用的函数句柄。
            % 它会自动检测模型是否包含 u 和 p，并相应地调整函数的输入参数。
            
            if nargin < 2
                func_name = 'f';
            end
            
            import casadi.*
            
            % 动态地构建输入参数列表和名称列表
            inputs = {obj.t_sym, obj.x_sym};
            input_names = {'t', 'x'};
            
            if obj.nu > 0
                inputs{end+1} = obj.u_sym;
                input_names{end+1} = 'u';
            end
            if obj.np > 0
                inputs{end+1} = obj.p_sym;
                input_names{end+1} = 'p';
            end
            
            % 创建CasADi函数
            dyn_fun = Function(func_name, inputs, {obj.ode_expr}, input_names, {'ode'});
        end

        function g_dot_fun = get_lie_derivative_function(obj, g_expr, func_name)
            % GET_LIE_DERIVATIVE_FUNCTION - 计算表达式 g(x) 沿系统动力学的时间导数。
            %
            % 这个方法实现了李导数的计算:
            %   d(g)/dt = (∂g/∂x) * f(t, x, u, p)
            %
            % 这在处理高阶路径约束或实现某些非线性控制策略时非常有用。
            %
            % 输入:
            %   g_expr - (casadi.SX) 一个主要依赖于状态 x 的CasADi符号表达式。
            %   func_name - (char, optional) 返回的CasADi函数的名称。
            %
            % 输出:
            %   g_dot_fun - (casadi.Function) 一个可调用的函数，其签名与
            %               get_dynamics_function 返回的函数一致，用于计算 d(g)/dt。

            if nargin < 3
                func_name = 'g_dot';
            end

            import casadi.*

            % --- 阶段一: 符号计算 ---
            % 验证 g_expr 是否为 casadi.SX 类型
            if ~isa(g_expr, 'casadi.SX')
                error('输入表达式 g_expr 必须是 casadi.SX 类型。');
            end

            % 计算 g(x) 关于 x 的雅可比矩阵 (梯度)
            J_g = jacobian(g_expr, obj.x_sym);
            
            % 获取系统动力学 f(t,x,u,p)
            f_expr = obj.ode_expr;
            
            % 计算李导数表达式: (∂g/∂x) * f
            lie_derivative_expr = J_g * f_expr + jacobian(g_expr, obj.t_sym);

            % --- 阶段二: 编译为 CasADi 函数 ---
            % 复用 get_dynamics_function 的逻辑来确保 API 一致性
            inputs = {obj.t_sym, obj.x_sym};
            input_names = {'t', 'x'};
            
            if obj.nu > 0
                inputs{end+1} = obj.u_sym;
                input_names{end+1} = 'u';
            end
            if obj.np > 0
                inputs{end+1} = obj.p_sym;
                input_names{end+1} = 'p';
            end
            
            % 创建与 get_dynamics_function 签名相同的 CasADi 函数
            g_dot_fun = Function(func_name, inputs, {lie_derivative_expr}, input_names, {'g_dot'});
        end
    end
end
