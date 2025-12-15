function [f, gradf] = objective_wrapper(w_vec, obj_fcn, grad_obj_fcn)
    % 包装目标函数

    % 计算目标函数值
    f = full(obj_fcn(w_vec));
    
    % 如果 fmincon 请求，则计算梯度
    if nargout > 1
        gradf = full(grad_obj_fcn(w_vec));
    end
end