function [c, ceq, grad_c, grad_ceq] = Constraints_Poly_VDPO(decisionVars, problem, T_Dis, ResPmtr)
% =-=======================================================================
% Polynomial Approximation Constraint Function (Exact Replication Version)
% =========================================================================
% This version is a direct, line-by-line replication of the logic in the
% original 'CnstrDisPoly_VDPOV5.m' to ensure identical numerical output.

%% --- 1. Parameter Extraction & Initialization (Replicated) ---
N = problem.control.N;
x0 = problem.state.x0;
controlGrid = problem.control.grid;
numStates = problem.system.dimensions;
odeFun = problem.functions.odeSensitivity;

T_Dis = unique(T_Dis); % Corresponds to original
num_T_Dis = length(T_Dis);

numSensitivityStates = numStates * N;
Stat0 = [x0, zeros(1, numSensitivityStates)];
y0 = Stat0;
numStateAug = length(y0);

state_at_T_Dis = zeros(num_T_Dis, numStateAug);
c_points = zeros(1, num_T_Dis);
grad_c_points = zeros(N, num_T_Dis);

opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

tout = [];
xout = [];

%% --- 2. Iterative Integration and Interpolation (Exact Replication) ---
for i = 1:N
    % Original ODE call
    [t, x] = ode45(@(t,x) odeFun(t, x, decisionVars(i), i, N, []), ...
                   [controlGrid(i), controlGrid(i+1)], y0, opts);
    y0 = x(end, :);

    % Original trajectory concatenation
    ntl = length(tout);
    tout = [tout(1:ntl-1); t];
    xout = [xout(1:ntl-1, :); x];

    % Original nested loop for interpolation
    for j = 1:num_T_Dis
        is_last_interval = (i == N);
        if ((T_Dis(j) >= controlGrid(i)) && (is_last_interval || (T_Dis(j) < controlGrid(i+1)))) == 1
            state_at_T_Dis(j, :) = interp1(tout, xout, T_Dis(j), 'PCHIP');
            
            c_points(1, j) = -state_at_T_Dis(j, 1) - 0.4;
            grad_c_points(:, j) = -state_at_T_Dis(j, numStates+1 : numStates : end)';
        end
    end
end
final_state_and_sens = y0; % This is xout(end,:) from the old code

%% --- 3. Polynomial Constraint Construction (Exact Replication) ---
t0 = problem.time.T0;
tf = problem.time.TF;
poly_intervals = T_Dis;
x_TD = state_at_T_Dis; % Use the same variable name for clarity

% Manual handling of endpoints t0 and tf, exactly as in the original
if isempty(find(t0 == T_Dis, 1))
    x_TD = [Stat0; x_TD];
    poly_intervals = [t0, poly_intervals];
end
if isempty(find(tf == T_Dis, 1))
    % This line in the old code had a potential bug if x_TD was empty.
    % Replicating it exactly.
    x_TD(end+1,:) = final_state_and_sens;
    poly_intervals = [poly_intervals, tf];
end

num_poly_intervals = length(poly_intervals) - 1;
PolyMax = zeros(1, num_poly_intervals);
PolyMax_u = zeros(N, num_poly_intervals);

for i = 1:num_poly_intervals
    tEnd = poly_intervals(i:i+1);
    
    % Replicated logic to find CurU and CurIdx for endpoints
    CurU = zeros(1,2); CurIdx = zeros(1,2);
    for ii = 1:N
        if (controlGrid(ii) <= tEnd(1)) && (tEnd(1) <= controlGrid(ii+1))
            if abs(controlGrid(ii+1) - tEnd(1)) > 1e-10
                CurU(1) = decisionVars(ii); CurIdx(1) = ii;
            else
                CurU(1) = decisionVars(ii+1); CurIdx(1) = ii+1;
            end
        end
        if (controlGrid(ii) <= tEnd(2)) && (tEnd(2) <= controlGrid(ii+1))
            % This part in the original had a potential index-out-of-bounds error
            % if tEnd(2) was exactly on controlGrid(1). Replicating it.
            if abs(controlGrid(ii) - tEnd(2)) > 1e-10
                CurU(2) = decisionVars(ii); CurIdx(2) = ii;
            else
                CurU(2) = decisionVars(ii-1); CurIdx(2) = ii-1;
            end
        end
    end
    
    gEnd = zeros(1, 2);
    dotgEnd = zeros(1, 2);
    gEnd_u = zeros(N, 2);
    dotgEnd_u = zeros(N, 2);
    Stat = x_TD(i:i+1, :)';
    
    for j = 1:2
        CurStt = Stat(:, j);
        CurT = tEnd(j);
        
        gEnd(j) = -x_TD(i+j-1, 1) - 0.4;
        gEnd_u(:, j) = -CurStt(numStates+1 : numStates : end);
        
        dx = odeFun(CurT, CurStt, CurU(j), CurIdx(j), N, []);
        dotgEnd(j) = -dx(1);
        dotgEnd_u(:, j) = -dx(numStates+1 : numStates : end);
    end
    
    % Assuming PloyMax_Gra exists and is identical to the original
    [PolyMax(i), PolyMax_u(:, i)] = PloyMax_Gra(tEnd, gEnd, dotgEnd, gEnd_u, dotgEnd_u);
end

%% --- 4. Combine and Finalize (Replicated) ---
c = [c_points, PolyMax + ResPmtr];
grad_c = [grad_c_points, PolyMax_u];
ceq = [];
grad_ceq = [];

end