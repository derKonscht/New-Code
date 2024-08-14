% Main loop
for t = 1:T
    % Update production functions
    Y_N(t) = (Y_cN^((epsilon-1)/epsilon) + Y_dN^((epsilon-1)/epsilon))^(epsilon/(epsilon-1));
    Y_cN = L_cN^(1-alpha_N) * integral(A_cN * x_cN^(1-alpha_N), 0, 1);
    Y_dN = R_N^alpha_2 * L_dN^(1-alpha_N) * integral(A_dN * x_dN^(1-alpha_N), 0, 1);
    Y_S(t) = (Y_cS^((epsilon-1)/epsilon) + Y_dS^((epsilon-1)/epsilon))^(epsilon/(epsilon-1));
    Y_cS = L_cS^(1-alpha_S) * integral(A_cS * x_cS^(1-alpha_S), 0, 1);
    Y_dS = R_S^alpha_2 * L_dS^(1-alpha_S) * integral(A_dS * x_dS^(1-alpha_S), 0, 1);

    % Update environmental quality
    S_N(t+1) = min(max(S_N(t) - qsi * Y_dN(t) + delta * S_N(t), 0), S_bar);
    S_S(t+1) = min(max(S_S(t) - qsi * Y_dS(t) + delta * S_S(t), 0), S_bar);

    % Update innovation
    A_cN(t+1) = (1 + eta_cN * gamma) * A_cN(t);
    A_dN(t+1) = (1 + eta_dN * gamma) * A_dN(t);
    A_cS(t+1) = kappa_c * A_cN(t) + (1 - kappa_c) * A_cS(t);
    A_dS(t+1) = kappa_d * A_dN(t) + (1 - kappa_d) * A_dS(t);

    % Calculate consumption
    C_N(t) = Y_N(t) - psi * (integral(x_cN) + integral(x_dN)) - c(Q_N) * R_N;
    C_S(t) = Y_S(t) - psi * (integral(x_cS) + integral(x_dS)) - c(Q_S) * R_S;

    % Utility and welfare
    U_N = @(C_N, S_N) u(C_N, S_N);
    U_S = @(C_S, S_S) u(C_S, S_S);
    W_N = sum((1/(1+rho)).^t .* U_N(C_N(t), S_N(t)));
    W_S = sum((1/(1+rho)).^t .* U_S(C_S(t), S_S(t)));
end

