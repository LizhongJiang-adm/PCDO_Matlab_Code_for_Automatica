function [c, ceq, gradc, gradceq] = nonlcon_wrapper(w_vec, nonlcon_fcn, grad_nonlcon_fcn, shift_value)
    % 包装非线性约束函数
    
    % 计算约束值
    % To ensure the final solution strictly satisfies g(u) <= 0 rather than g(u) <= Tol,
    % we shift the constraints passed to fmincon by its own tolerance.
    [c_cas, ceq_cas] = nonlcon_fcn(w_vec);
    c = full(c_cas) +  shift_value;
    ceq = full(ceq_cas);
    
    % 如果 fmincon 请求，则计算梯度
    if nargout > 2
        [gradc_cas, gradceq_cas] = grad_nonlcon_fcn(w_vec);
        gradc = full(gradc_cas);
        gradceq = full(gradceq_cas);
    end
end