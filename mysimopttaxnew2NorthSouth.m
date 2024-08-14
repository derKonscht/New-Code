function U = mysimopttaxnew2NorthSouth(x, Ac0_north, Ad0_north, S0_north, Ac0_south, Ad0_south, S0_south, Ac_new, Ad_new, Ac_old, Ad_old, kappa)
    global rho sigma psi alpha gamma eta_d_north eta_c_north eta_d_south eta_c_south qsi epsilon delta numsim S_bar

    % Separate vectors for North and South
    s_c_north = reshape(x(1:numsim), [1, numsim]);
    s_c_south = reshape(x(numsim+1:2*numsim), [1, numsim]);
    tau_north = reshape(x(2*numsim+1:3*numsim), [1, numsim]);
    tau_south = reshape(x(3*numsim+1:end), [1, numsim]);

    % Initialize vectors
    A_c_north = zeros(numsim, 1);
    A_d_north = zeros(numsim, 1);
    C_north = zeros(numsim, 1);
    S_north = zeros(numsim, 1);
    A_c_south = zeros(numsim, 1);
    A_d_south = zeros(numsim, 1);
    C_south = zeros(numsim, 1);
    S_south = zeros(numsim, 1);

    % Initial values
    S_north(1) = S0_north;
    S_south(1) = S0_south;
    A_c_north(1) = (1 + gamma * eta_c_north * s_c_north(1)) * Ac0_north;
    A_d_north(1) = (1 + gamma * eta_d_north * (1 - s_c_north(1))) * Ad0_north;
    A_c_south(1) = kappa * Ac_new + (1 - kappa) * Ac_old;
    A_d_south(1) = kappa * Ad_new + (1 - kappa) * Ad_old;

    % Calculate consumption for both North and South
    C_north(1) = (alpha / psi)^(alpha / (1 - alpha)) * A_c_north(1) * A_d_north(1) / ((1 + tau_north(1))^(1 - epsilon) * A_c_north(1)^((1 - alpha) * (1 - epsilon)) + A_d_north(1)^((1 - alpha) * (1 - epsilon)))^(1 / ((1 - alpha) * (1 - epsilon))) * (1 - alpha + tau_north(1) * A_c_north(1)^((1 - alpha) * (1 - epsilon)) / (A_c_north(1)^((1 - alpha) * (1 - epsilon)) + (1 + tau_north(1))^epsilon * A_d_north(1)^((1 - alpha) * (1 - epsilon))));
    C_south(1) = (alpha / psi)^(alpha / (1 - alpha)) * A_c_south(1) * A_d_south(1) / ((1 + tau_south(1))^(1 - epsilon) * A_c_south(1)^((1 - alpha) * (1 - epsilon)) + A_d_south(1)^((1 - alpha) * (1 - epsilon)))^(1 / ((1 - alpha) * (1 - epsilon))) * (1 - alpha + tau_south(1) * A_c_south(1)^((1 - alpha) * (1 - epsilon)) / (A_c_south(1)^((1 - alpha) * (1 - epsilon)) + (1 + tau_south(1))^epsilon * A_d_south(1)^((1 - alpha) * (1 - epsilon))));

    % Simulation
    for n = 2:numsim
        % Update North
        A_c_north(n) = (1 + gamma * eta_c_north * s_c_north(n)) * A_c_north(n - 1);
        A_d_north(n) = (1 + gamma * eta_d_north * (1 - s_c_north(n))) * A_d_north(n - 1);
        C_north(n) = (alpha / psi)^(alpha / (1 - alpha)) * A_c_north(n) * A_d_north(n) / ((1 + tau_north(n))^(1 - epsilon) * A_c_north(n)^((1 - alpha) * (1 - epsilon)) + A_d_north(n)^((1 - alpha) * (1 - epsilon)))^(1 / ((1 - alpha) * (1 - epsilon))) * (1 - alpha + tau_north(n) * A_c_north(n)^((1 - alpha) * (1 - epsilon)) / (A_c_north(n)^((1 - alpha) * (1 - epsilon)) + (1 + tau_north(n))^epsilon * A_d_north(n)^((1 - alpha) * (1 - epsilon))));
        
        % Update South
        A_c_south(n) = kappa * (1 + gamma * eta_c_north * s_c_north(n)) * Ac_new + (1 - kappa) * (1 + gamma * eta_c_north * s_c_north(n)) * Ac_old;
       
        % Initialize Y_d for both regions
Y_d_N = zeros(numsim,1);
Y_d_S = zeros(numsim,1);

for n = 2:numsim
    % Existing calculations for A_c, A_d, Q, C, etc.

    % Update Y_d_N and Y_d_S based on respective country parameters
    Y_d_N(n) = % (compute Y_d for North)
    Y_d_S(n) = % (compute Y_d for South)

    % Modify the environmental quality update to include both regions
    S(n) = min(max(0.00000000000000001, -qsi * (Y_d_N(n-1) + Y_d_S(n-1)) + (1 + delta) * S(n-1)), S_bar);
end

