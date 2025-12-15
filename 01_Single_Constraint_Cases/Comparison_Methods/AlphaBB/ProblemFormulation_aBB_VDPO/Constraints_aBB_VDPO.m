function [c, ceq, grad_c, grad_ceq] = Constraints_aBB_VDPO(decisionVars, problem, eta)
% =========================================================================
% aBB Method Constraint Function for VDPO (Standardized)
% =========================================================================
% Description:
% Computes the nonlinear inequality (c) and equality (ceq) constraints
% for the aBB method. This includes the凹化(convexified) path constraint
% and smoothed KKT complementarity conditions. This version is self-contained.
%
% Inputs:
%   decisionVars - The composite vector of decision variables [u; t_k; gamma_L; gamma_U].
%   problem      - The problem definition struct.
%
% Outputs:
%   c, ceq         - Inequality and equality constraint vectors.
%   grad_c, grad_ceq - Gradients of the constraints w.r.t. all decision variables.
% =========================================================================

%% --- 1. Extract Parameters from Problem Struct ---
x0 = problem.state.x0;
N = problem.control.N;
T = problem.time.TF;
controlGrid = problem.control.grid;
numStates = problem.system.dimensions;
odeFun = problem.functions.odeSensitivity;

% aBB specific parameters
M_eta = length(eta) - 1;
myalpha = problem.algorithm_params.myalpha; % Convexification parameter
tau = problem.algorithm_params.tau;

% Handle to the smoothed complementarity function
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
numSensitivityStates = numStates * N;
y0 = [x0, zeros(1, numSensitivityStates)];

% Initialize arrays to store the full concatenated trajectory
tout = [];
xout = [];
opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9); % Use high precision for accuracy

for i = 1:N
    t_start = controlGrid(i);
    t_end = controlGrid(i+1);
    time_span = [t_start, t_end];
    u_i = u_control(i); % Control for this interval
    
    % The ODE handle for this specific interval
    ode_handle_interval = @(t, x) odeFun(t, x, u_i, i, N, [t_start, t_end]);
    
    % Integrate over the small interval
    [t_interval, x_interval] = ode45(ode_handle_interval, time_span, y0, opts);
    
    % Update initial state for the next interval
    y0 = x_interval(end, :);
    
    % Append results, avoiding duplicate time points
    if isempty(tout)
        tout = t_interval;
        xout = x_interval;
    else
        tout = [tout; t_interval(2:end)];
        xout = [xout; x_interval(2:end, :)];
    end
end

% Interpolate to get state and sensitivity values at each t_k
state_and_sens_at_tk = interp1(tout, xout, t_k, 'PCHIP');

%% --- 5. Build Constraints and Gradients for Each Interval ---
for i = 1:M_eta
    % Find the corresponding control u_star for the interval [eta(i), eta(i+1)]
    control_idx = find(controlGrid <= eta(i) + 1e-10, 1, 'last');
    u_star = u_control(control_idx);
    
    % Extract state values at t_k(i)
    x_at_tk = state_and_sens_at_tk(i, 1:numStates);
    x1 = x_at_tk(1); x2 = x_at_tk(2);

    % --- A. Inequality Constraint c(i) ---
    g_tk = -x1 - 0.4;
    c(i) = g_tk + (myalpha/2) * (t_k(i) - eta(i)) * (eta(i+1) - t_k(i));
    
    % --- B. Equality Constraints ceq(i), ceq(M+i), ceq(2M+i) ---
    g_dot_tk = -(1 - x2^2) * x1 + x2 - u_star;
    ceq(i) = g_dot_tk + myalpha * ((eta(i)+eta(i+1))/2 - t_k(i)) + gamma_L(i) - gamma_U(i);
    ceq(M_eta+i) = phi_fun(gamma_L(i), t_k(i) - eta(i));
    ceq(2*M_eta+i) = phi_fun(gamma_U(i), eta(i+1) - t_k(i));
    
    % --- C. Compute Gradients ---
    % C.1 Gradient w.r.t. control variables u_p (p=1...N)
    for p = 1:N
        sens_at_tk = state_and_sens_at_tk(i, numStates*p+1 : numStates*(p+1));
        x1_q_p = sens_at_tk(1); x2_q_p = sens_at_tk(2);
        
        % Grad of c(i)
        grad_c(p, i) = -x1_q_p;
        
        % Grad of ceq(i)
        dg_dot_dx1 = x2^2 - 1.0;
        dg_dot_dx2 = 2*x1*x2 + 1.0;
        grad_ceq(p, i) = dg_dot_dx1 * x1_q_p + dg_dot_dx2 * x2_q_p;
        if p == control_idx
            grad_ceq(p, i) = grad_ceq(p, i) - 1; % Direct derivative
        end
    end
    
    % C.2 Gradient w.r.t. time variables t_k(j) (j=1...M_eta)
    % Gradient is non-zero only for j=i (diagonal structure)
    % Note: Derivative of interp1 is complex. Here we assume chain rule applies directly to g_dot_tk
    grad_c(N+i, i) = g_dot_tk + myalpha * ((eta(i)+eta(i+1))/2 - t_k(i));
    
    % g_dot_dot_tk = g_dot_dot_fun(x1, x2, x_at_tk(3), u_star, t_k(i));
    g_dot_dot_tk = (1 - x2^2) * g_dot_tk + 2*x1^2*x2 + x1; % Note: -(x2^2-1) = (1-x2^2)
    grad_ceq(N+i, i) = g_dot_dot_tk - myalpha;
    
    sqrt_term_L = sqrt(gamma_L(i)^2 + (t_k(i)-eta(i))^2 + 2*tau^2);
    grad_ceq(N+i, M_eta+i) = 1 - (t_k(i)-eta(i)) / sqrt_term_L;
    
    sqrt_term_U = sqrt(gamma_U(i)^2 + (eta(i+1)-t_k(i))^2 + 2*tau^2);
    grad_ceq(N+i, 2*M_eta+i) = - (1 - (eta(i+1)-t_k(i)) / sqrt_term_U);

    % C.3 Gradient w.r.t. gamma_L(j)
    grad_ceq(N+M_eta+i, i) = 1;
    grad_ceq(N+M_eta+i, M_eta+i) = 1 - gamma_L(i) / sqrt_term_L;
    
    % C.4 Gradient w.r.t. gamma_U(j)
    grad_ceq(N+2*M_eta+i, i) = -1;
    grad_ceq(N+2*M_eta+i, 2*M_eta+i) = 1 - gamma_U(i) / sqrt_term_U;
end

end