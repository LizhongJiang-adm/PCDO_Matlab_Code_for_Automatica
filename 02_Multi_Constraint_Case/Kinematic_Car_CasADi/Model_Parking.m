% Model_Parking.m
classdef Model_Parking < DynamicModel
    % Model_Parking - 泊车问题的自行车运动学模型。
    
    properties (Constant)
        % vehicle_params - 车辆的几何参数
        % 注意: 根据要求，n 为前悬，m 为后悬
        vehicle_params = struct(...
            'l', 2.8, ...   % 轴距
            'n', 0.96, ...  % 前悬 (Front Overhang)
            'm', 0.929, ... % 后悬 (Rear Overhang)
            'b', 0.971 ...  % 半车宽
        );
    end
    
    methods
        function obj = Model_Parking()
            obj@DynamicModel();
        end
        
        function setup_model(obj)
            import casadi.*
            
            % x_sym (5x1): [px, py, theta, v, phi]
            obj.x_sym = SX.sym('x_state', 5, 1);
            
            % u_sym (2x1): [a, omega]
            obj.u_sym = SX.sym('u_signal', 2, 1);
            
            x = obj.x_sym; u = obj.u_sym;
            px=x(1); py=x(2); theta=x(3); v=x(4); phi=x(5);
            a=u(1); omega=u(2);
            l = obj.vehicle_params.l;
            
            obj.ode_expr = [v * cos(theta);
                            v * sin(theta);
                            (v / l) * tan(phi);
                            a;
                            omega];
        end
        
    end
end