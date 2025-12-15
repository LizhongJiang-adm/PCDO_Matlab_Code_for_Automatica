function [obj_fcn, grad_obj_fcn, nonlcon_fcn, grad_nonlcon_fcn, ...
    initial_guess, lbw, ubw] = NmlzNLP(nlp)

        w_sym = nlp.x;
        Idxeq = find(nlp.lbg==nlp.ubg);
        ceq_sym = nlp.g(Idxeq) - nlp.lbg(Idxeq);
        IdxEffIeqL = find(nlp.lbg ~= nlp.ubg   &   nlp.lbg ~= -Inf);
        IdxEffIeqU = find(nlp.lbg ~= nlp.ubg   &   nlp.ubg ~= Inf);
        cieq_sym = [ - nlp.g(IdxEffIeqL) + nlp.lbg(IdxEffIeqL) ; nlp.g(IdxEffIeqU) - nlp.ubg(IdxEffIeqU)];
        
        jac_c_sym = jacobian(cieq_sym, w_sym)';
        jac_ceq_sym = jacobian(ceq_sym, w_sym)';
        
        nonlcon_fcn = casadi.Function('nonlcon_fcn', {w_sym}, {cieq_sym, ceq_sym});
        grad_nonlcon_fcn = casadi.Function('grad_nonlcon_fcn', {w_sym}, {jac_c_sym, jac_ceq_sym});
        
        J = nlp.f;
        grad_J_sym = jacobian(J, w_sym)'; % fmincon 需要列向量梯度
        
        obj_fcn = casadi.Function('obj_fcn', {w_sym}, {J});
        grad_obj_fcn = casadi.Function('grad_obj_fcn', {w_sym}, {grad_J_sym});
        lbw = nlp.lbx;  ubw = nlp.ubx;
    
        initial_guess = nlp.x0;

end