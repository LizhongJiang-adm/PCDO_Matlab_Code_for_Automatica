function result = phi_tau_refactored(a, b, problem)
% =========================================================================
% Smoothed Complementarity Function for aBB Method (Refactored)
% =========================================================================
    % This function approximates the complementarity condition a*b=0.
    % It uses a parameter tau, which should be defined in the problem struct.
    tau = problem.algorithm_params.tau;
    result = a + b - sqrt(a.^2 + b.^2 + 2*tau^2);
end