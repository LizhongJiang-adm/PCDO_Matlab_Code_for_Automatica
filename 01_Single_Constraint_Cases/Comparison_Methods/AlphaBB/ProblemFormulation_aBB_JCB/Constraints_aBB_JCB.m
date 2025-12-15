function [c, ceq, grad_c, grad_ceq] = Constraints_aBB_JCB(decisionVars, problem, eta)
% =========================================================================
% aBB Method Constraint Function for JCB Problem (Standardized)
% =========================================================================
% ... (Description, Inputs, Outputs are the same as VDPO version)

%% --- 1. Extract Parameters from Problem Struct ---
x0 = problem.state.x0;
N = problem.control.N;
T = problem.time.TF;
controlGrid = problem.control.grid;
numStates = problem.system.dimensions;
odeFun = problem.functions.odeSensitivity;

% aBB specific parameters
M_eta = length(eta) - 1;
myalpha = problem.algorithm_params.myalpha;
tau = problem.algorithm_params.tau;
phi_fun = @(a,b) phi_tau_refactored(a,b,problem);

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
% This section is identical to the VDPO version
numSensitivityStates = numStates * N;
y0 = [x0, zeros(1, numSensitivityStates)];
tout = []; xout = [];
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

for i = 1:N
    t_start = controlGrid(i);
    t_end = controlGrid(i+1);
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
    x2 = x_at_tk(2);

    % --- A. Inequality Constraint c(i) for JCB PROBLEM ---
    g_tk = x2 - 8.0 * (t_k(i) - 0.5)^2 + 0.5;
    c(i) = g_tk + (myalpha/2) * (t_k(i) - eta(i)) * (eta(i+1) - t_k(i));
    
    % --- B. Equality Constraints ceq(i) for JCB PROBLEM ---
    g_dot_tk = -16.0 * t_k(i) + u_star - x2 + 8.0;
    ceq(i) = g_dot_tk + myalpha * ((eta(i)+eta(i+1))/2 - t_k(i)) + gamma_L(i) - gamma_U(i);
    ceq(M_eta+i) = phi_fun(gamma_L(i), t_k(i) - eta(i));
    ceq(2*M_eta+i) = phi_fun(gamma_U(i), eta(i+1) - t_k(i));
    
    % --- C. Compute Gradients for JCB PROBLEM ---
    % C.1 Gradient w.r.t. control variables u_p (p=1...N)
    for p = 1:N
        % For JCB, numStates=2 or 3 depending on cost. Assuming 3 for generality
        sens_at_tk = state_and_sens_at_tk(i, numStates*p+1 : numStates*(p+1));
        x2_q_p = sens_at_tk(2);
        
        grad_c(p, i) = x2_q_p; % Grad of g_tk w.r.t u_p is dx2/du_p
        
        grad_ceq(p, i) = -x2_q_p; % Grad of g_dot_tk w.r.t u_p
        if p == control_idx
            grad_ceq(p, i) = grad_ceq(p, i) + 1; % Direct derivative
        end
    end
    
    % C.2 Gradient w.r.t. time variables t_k(j) (j=1...M_eta)
    % grad_c(N+i, i) = d/dt_k(g_tk) + d/dt_k(凹化项)
    grad_c(N+i, i) = g_dot_tk + myalpha * ((eta(i)+eta(i+1))/2 - t_k(i)); % Corrected version
    
    g_dot_dot_tk = -u_star + x2 - 16.0; % Using formula from original JCB code
    grad_ceq(N+i, i) = g_dot_dot_tk - myalpha;
    
    sqrt_term_L = sqrt(gamma_L(i)^2 + (t_k(i)-eta(i))^2 + 2*tau^2);
    grad_ceq(N+i, M_eta+i) = 1 - (t_k(i)-eta(i)) / sqrt_term_L;
    
    sqrt_term_U = sqrt(gamma_U(i)^2 + (eta(i+1)-t_k(i))^2 + 2*tau^2);
    grad_ceq(N+i, 2*M_eta+i) = - (1 - (eta(i+1)-t_k(i)) / sqrt_term_U);

    % C.3 & C.4: Gradients w.r.t. gamma_L and gamma_U (these are problem-independent)
    grad_ceq(N+M_eta+i, i) = 1;
    grad_ceq(N+M_eta+i, M_eta+i) = 1 - gamma_L(i) / sqrt_term_L;
    
    grad_ceq(N+2*M_eta+i, i) = -1;
    grad_ceq(N+2*M_eta+i, 2*M_eta+i) = 1 - gamma_U(i) / sqrt_term_U;
end

end