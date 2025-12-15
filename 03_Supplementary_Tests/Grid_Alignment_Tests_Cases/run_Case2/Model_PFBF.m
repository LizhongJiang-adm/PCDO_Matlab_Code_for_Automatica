
classdef  Model_PFBF < DynamicModel
    
    properties

    end
    
    methods
        
        function obj = Model_PFBF()
            obj@DynamicModel(); 

        end
        
        function setup_model(obj)
            import casadi.*
            
            % --- 1. 定义符号化的模型参数 (p_sym) ---
            % 只包含车辆的几何参数。

            % --- 2. 定义状态、控制和时间变量 ---
            obj.x_sym = SX.sym('x_state', 4, 1); % [px, py, v, theta, phi]
            obj.u_sym = SX.sym('u_signal', 1, 1); % [a, omega]
            obj.t_sym = SX.sym('t_time');       % 真实时间 t
            
            % --- 3. 构建核心ODE表达式 (dx/dt) ---
            x = obj.x_sym;
            u = obj.u_sym;
            % p = obj.p_sym;
            
            kL=0.006;mu=0.11;
            Yxs=0.47;
            theta=0.004;Yp=1.2;kI=0.1;Mx=0.029;Sl=400;kXP=0.01;
            kP=0.0001;


            obj.ode_expr = [(mu*x(1)*x(2))/(x(2) + kL*x(1)) - (u(1)*x(1))/x(4);
                            (u(1)*(Sl - x(2)))/x(4) - Mx*x(1) - (theta*x(1)*x(2))/(Yp*(kP + x(2) + x(2)^2/kI)) - (mu*x(1)*x(2))/(Yxs*(x(2) + kL*x(1)));
                            (theta*x(1)*x(2))/(kP + x(2) + x(2)^2/kI) - (u(1)*x(3))/x(4) - kXP*x(3);
                            u(1);];

        end
        
    end
end