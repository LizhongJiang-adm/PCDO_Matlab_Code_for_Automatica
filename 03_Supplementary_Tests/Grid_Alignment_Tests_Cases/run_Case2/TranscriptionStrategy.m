% TranscriptionStrategy.m
classdef (Abstract) TranscriptionStrategy < handle
    % TranscriptionStrategy - 从连续OCP到离散NLP的转译策略的抽象基类。
    %
    % 这是一个契约，所有具体的离散化方法(如多重打靶法、直接配置法)
    % 都应继承此类，并实现其抽象方法。
    %
    % 它的核心职责是接收一个 OptimalControlProblem 对象和一个时间网格，
    % 并将其“翻译”成一个结构化的、可被NLP求解器理解的离散问题。
    
    properties (SetAccess = protected)
        % --- 输入 ---
        ocp                 % 存储传入的 OptimalControlProblem 对象
        grid_points         % 存储离散化的时间网格点
        
        % --- NLP 输出 (由子类在 transcribe 方法中填充) ---
        % 这些属性构成了与具体离散化结构无关的、扁平化的NLP问题定义
        
        nlp_vars_sym        % 所有决策变量“压平”后的一个CasADi符号向量
        nlp_vars_lb         % 决策变量向量的数值下界
        nlp_vars_ub         % 决策变量向量的数值上界
        nlp_initial_guess   % 决策变量向量的数值初始猜测
        
        nlp_objective       % 最终的标量目标函数表达式
        nlp_constraints     % 所有约束拼接成的一个CasADi向量表达式
        nlp_constraints_lb  % 约束向量的数值下界
        nlp_constraints_ub  % 约束向量的数值上界
    end
    
    methods (Abstract)
        % --- 抽象方法：子类的强制性任务 ---
        
        % transcribe - 执行从连续到离散的转译过程
        %
        % 这是所有具体策略类都必须实现的“核心引擎”。它负责:
        % 1. 根据具体策略(如配置法)创建所有离散化的符号决策变量。
        % 2. 将所有这些变量“压平”成单一向量，并存入 obj.nlp_vars_sym。
        % 3. 根据OCP中的边界和策略结构，创建并填充扁平化的数值边界向量
        %    obj.nlp_vars_lb 和 obj.nlp_vars_ub。
        % 4. 根据策略结构，创建并填充扁平化的初始猜测向量 obj.nlp_initial_guess。
        % 5. 构建离散化的目标函数，并存入 obj.nlp_objective。
        % 6. 构建所有离散化的约束(动力学、路径、边界等)，拼接成单一向量，
        %    并与它们的数值边界一同存入 obj.nlp_constraints, 
        %    obj.nlp_constraints_lb, obj.nlp_constraints_ub。
        transcribe(obj);
    end
    
    methods
        function obj = TranscriptionStrategy(ocp, grid_points)
            % 构造函数
            % 接收一个OCP对象和一个灵活的时间网格。
            
            if ~isa(ocp, 'OptimalControlProblem')
                error('输入必须是 OptimalControlProblem 类的对象。');
            end
            if ~isnumeric(grid_points) || ~isvector(grid_points) || any(diff(grid_points) <= 0)
                error('grid_points 必须是一个单调递增的数值向量。');
            end
            
            obj.ocp = ocp;
            obj.grid_points = grid_points;
        end
        
        function nlp_problem = get_nlp_problem(obj)
            % get_nlp_problem - 提供一个标准接口来获取最终的NLP问题
            %
            % 这个方法将转译后的组件打包成一个可被CasADi求解器
            % (如 nlpsol) 直接使用的标准结构体。
            
            if isempty(obj.nlp_objective) || isempty(obj.nlp_constraints)
                error('OCP尚未被转译。请在调用此方法前先调用 transcribe() 方法。');
            end
            
            nlp_problem = struct(...
                'x',   obj.nlp_vars_sym, ...      % (x) 决策变量符号向量
                'f',   obj.nlp_objective, ...     % (f) 目标函数表达式
                'g',   obj.nlp_constraints, ...   % (g) 约束函数表达式向量
                'lbx', obj.nlp_vars_lb, ...       % (lbx) 决策变量下界
                'ubx', obj.nlp_vars_ub, ...       % (ubx) 决策变量上界
                'lbg', obj.nlp_constraints_lb, ...% (lbg) 约束下界
                'ubg', obj.nlp_constraints_ub, ...% (ubg) 约束上界
                'x0',  obj.nlp_initial_guess ...  % (x0) 决策变量初始猜测
            );
        end
    end
end