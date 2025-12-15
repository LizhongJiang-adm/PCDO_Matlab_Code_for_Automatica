function [c, ceq, grad_c, grad_ceq] = Constraints_aBB_PFBF(decisionVars, problem, eta)
% =========================================================================
% aBB Method Constraint Function for PFBF Problem (Standardized)
% =========================================================================
% ... (Description, Inputs, Outputs are the same as other aBB versions)

%% --- 1. Extract Parameters from Inputs ---
x0 = problem.state.x0;
N = problem.control.N;
controlGrid = problem.control.grid;
numStates = problem.system.dimensions;
odeFun = problem.functions.odeSensitivity;

% aBB specific parameters
M_eta = length(eta) - 1;
myalpha = problem.algorithm_params.myalpha;
tau = problem.algorithm_params.tau;
phi_fun = @(a,b) phi_tau_refactored(a,b,problem);

% Physical parameters for the PFBF model
% It's better to pass these via 'problem' struct, but for now we define them locally.
kL=0.006; mu=0.11; Yxs=0.47; theta=0.004; Yp=1.2;
kI=0.1; Mx=0.029; Sl=400; kXP=0.01; kP=0.0001;

%% --- 2. Decompose the Composite Decision Variable Vector ---
num_total_vars = N + 3*M_eta;
u_control = decisionVars(1:N);
t_k = decisionVars(N+1 : N+M_eta);
gamma_L = decisionVars(N+M_eta+1 : N+2*M_eta);
gamma_U = decisionVars(N+2*M_eta+1 : N+3*M_eta);

%% --- 3. Initialize Outputs ---
c = zeros(M_eta, 1);
ceq = zeros(3*M_eta, 1);
grad_c = zeros(num_total_vars, M_eta);
grad_ceq = zeros(num_total_vars, 3*M_eta);

%% --- 4. Integrate System Dynamics and Sensitivities ---
% This section is identical to the standardized versions
numSensitivityStates = numStates * N;
y0 = [x0, zeros(1, numSensitivityStates)];
tout = []; xout = [];
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

for i = 1:N
    t_start = controlGrid(i); t_end = controlGrid(i+1);
    u_i = u_control(i);
    ode_handle_interval = @(t, x) odeFun(t, x, u_i, i, N, [t_start, t_end]);
    [t_interval, x_interval] = ode45(ode_handle_interval, [t_start, t_end], y0, opts);
    y0 = x_interval(end, :);
    if isempty(tout)
        tout = t_interval; xout = x_interval;
    else
        tout = [tout; t_interval(2:end)]; xout = [xout; x_interval(2:end, :)];
    end
end
state_and_sens_at_tk = interp1(tout, xout, t_k, 'PCHIP');

%% --- 5. Build Constraints and Gradients for Each Interval ---
for i = 1:M_eta
    control_idx = find(controlGrid <= eta(i) + 1e-10, 1, 'last');
    u_star = u_control(control_idx);
    
    x_at_tk = state_and_sens_at_tk(i, 1:numStates);
    x1 = x_at_tk(1); x2 = x_at_tk(2); x3 = x_at_tk(3); x4 = x_at_tk(4);

    % --- A. Inequality Constraint c(i) for PFBF PROBLEM ---
    g_tk = x2 - 0.5;
    c(i) = g_tk + (myalpha/2) * (t_k(i) - eta(i)) * (eta(i+1) - t_k(i));
    
    % --- B. Equality Constraints ceq(i) for PFBF PROBLEM ---
    term1 = (u_star .* (Sl - x2)) ./ x4;
    term2 = Mx .* x1;
    term3_num = theta .* x1 .* x2;
    term3_den = Yp .* (x2.^2./kI + x2 + kP);
    term4_num = mu .* x1 .* x2;
    term4_den = Yxs .* (x2 + kL .* x1);
    g_dot_tk = term1 - term2 - term3_num ./ term3_den - term4_num ./ term4_den;
    
    ceq(i) = g_dot_tk + myalpha * ((eta(i)+eta(i+1))/2 - t_k(i)) + gamma_L(i) - gamma_U(i);
    ceq(M_eta+i) = phi_fun(gamma_L(i), t_k(i) - eta(i));
    ceq(2*M_eta+i) = phi_fun(gamma_U(i), eta(i+1) - t_k(i));
    
    % --- C. Compute Gradients for PFBF PROBLEM ---
    % C.1 Gradient w.r.t. control variables u_p (p=1...N)
    for p = 1:N
        sens_at_tk_p = state_and_sens_at_tk(i, numStates*p+1 : numStates*(p+1));
        x1_q_p = sens_at_tk_p(1); x2_q_p = sens_at_tk_p(2);
        x3_q_p = sens_at_tk_p(3); x4_q_p = sens_at_tk_p(4);
        
        grad_c(p, i) = x2_q_p; % Grad of g_tk w.r.t u_p is dx2/du_p
        
        dg_dot_dx1 = (kL.*mu.*x1.*x2)./(Yxs.*(x2 + kL.*x1).^2) - (mu.*x2)./(Yxs.*(x2 + kL.*x1)) - (theta.*x2)./(Yp.*(x2.^2./kI + x2 + kP)) - Mx;
        dg_dot_dx2 = (mu.*x1.*x2)./(Yxs.*(x2 + kL.*x1).^2) - (mu.*x1)./(Yxs.*(x2 + kL.*x1)) - (theta.*x1)./(Yp.*(x2.^2./kI + x2 + kP)) - u_star./x4 + (theta.*x1.*x2.*((2.*x2)./kI + 1))./(Yp.*(x2.^2./kI + x2 + kP).^2);
        dg_dot_dx3 = 0;
        dg_dot_dx4 = -(u_star.*(Sl - x2))./x4.^2;
        
        grad_ceq(p, i) = dg_dot_dx1*x1_q_p + dg_dot_dx2*x2_q_p + dg_dot_dx3*x3_q_p + dg_dot_dx4*x4_q_p;
        
        if p == control_idx
            dg_dot_du = (Sl - x2)./x4;
            grad_ceq(p, i) = grad_ceq(p, i) + dg_dot_du;
        end
    end
    
    % C.2 Gradient w.r.t. time variables t_k(j) (j=1...M_eta)
    grad_c(N+i, i) = g_dot_tk + myalpha * ((eta(i)+eta(i+1))/2 - t_k(i));
    
    dot_u_star = 0; % Assuming piecewise constant control
    term_A = (u_star*x1)/x4 - (mu*x1*x2)/(x2 + kL*x1);
    term_B = Mx + (mu*x2)/(Yxs*(x2 + kL*x1)) + (theta*x2)/(Yp*(x2^2/kI + x2 + kP)) - (kL*mu*x1*x2)/(Yxs*(x2 + kL*x1)^2);
    term_C = Mx*x1 - (u_star*(Sl - x2))/x4 + (theta*x1*x2)/(Yp*(x2^2/kI + x2 + kP)) + (mu*x1*x2)/(Yxs*(x2 + kL*x1));
    term_D = u_star/x4 + (mu*x1)/(Yxs*(x2 + kL*x1)) + (theta*x1)/(Yp*(x2^2/kI + x2 + kP)) - (mu*x1*x2)/(Yxs*(x2 + kL*x1)^2) - (theta*x1*x2*((2*x2)/kI + 1))/(Yp*(x2^2/kI + x2 + kP)^2);
    term_E = (dot_u_star*(Sl - x2))/x4 - (u_star^2*(Sl - x2))/x4^2;
    g_dot_dot_tk = term_A*term_B + term_C*term_D + term_E;
    
    grad_ceq(N+i, i) = g_dot_dot_tk - myalpha;
    
    % The rest of the gradients for phi_tau are problem-independent and correct
    sqrt_term_L = sqrt(gamma_L(i)^2 + (t_k(i)-eta(i))^2 + 2*tau^2);
    grad_ceq(N+i, M_eta+i) = 1 - (t_k(i)-eta(i)) / sqrt_term_L;
    
    sqrt_term_U = sqrt(gamma_U(i)^2 + (eta(i+1)-t_k(i))^2 + 2*tau^2);
    grad_ceq(N+i, 2*M_eta+i) = - (1 - (eta(i+1)-t_k(i)) / sqrt_term_U);

    grad_ceq(N+M_eta+i, i) = 1;
    grad_ceq(N+M_eta+i, M_eta+i) = 1 - gamma_L(i) / sqrt_term_L;
    
    grad_ceq(N+2*M_eta+i, i) = -1;
    grad_ceq(N+2*M_eta+i, 2*M_eta+i) = 1 - gamma_U(i) / sqrt_term_U;
end

end