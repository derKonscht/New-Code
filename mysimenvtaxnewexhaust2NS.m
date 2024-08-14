function Resp = mysimenvtaxnewexhaust2NS(x, Ac0, Ad0, S0, Q0)
    global rho sigma psi alpha gamma eta_d_N eta_c_N eta_d_S eta_c_S qsi epsilon delta numsim alpha2 alpha1 k1 k2 kk phi1 phi S_bar 

    s_c = x(1:1:numsim);
    tau = x((numsim+1):1:(2*numsim));
    theta = x((2*numsim+1):1:(3*numsim));

    A_c_N = zeros(numsim,1);
    A_d_N = zeros(numsim,1);
    A_c_S = zeros(numsim,1);
    A_d_S = zeros(numsim,1);
    C_N = zeros(numsim,1);
    C_S = zeros(numsim,1);
    s_c_N = zeros(numsim,1);
    s_c_S = zeros(numsim,1); 
    S = zeros(numsim,1);
    q = zeros(numsim,1);
    S(1) = S0; 
    Q = zeros(numsim,1);
    Q(1) = Q0;
%%%Initial values
S(1) = S0; 
A_c_N(1) = (1 + gamma * eta_c_N * s_c_N(1)) * Ac0;
A_d_N(1) = (1 + gamma * eta_d_N * (1 - s_c_N(1))) * Ad0;
A_c_S(1) = (kk* s_c_S(1) * A_c0) + ((1- kk)* (1-s_c_S(1))* Ac0;
A_d_S(1) = (kk* s_c_S(1) * A_c0) + ((1- kk)* (1-s_c_S(1))* Ad0;
%Consumption in the base paper is B.6 in the Appendix, check and change if
%necessary
C_N(1) = (alpha / psi)^(alpha / (1 - alpha)) * A_c_N(1) * A_d_N(1) / ((1 + tau(1))^(1 - epsilon) * A_c_N(1)^((1 - alpha) * (1 - epsilon)) + A_d_N(1)^((1 - alpha) * (1 - epsilon)))^(1 / ((1 - alpha) * (1 - epsilon))) * (1 - alpha + tau(1) * A_c_N(1)^((1 - alpha) * (1 - epsilon)) / (A_c_N(1)^((1 - alpha) * (1 - epsilon)) + (1 + tau(1))^epsilon * A_d_N(1)^((1 - alpha) * (1 - epsilon))));
C_S(1) = (alpha / psi)^(alpha / (1 - alpha)) * A_c_S(1) * A_d_S(1) / ((1 + tau(1))^(1 - epsilon) * A_c_S(1)^((1 - alpha) * (1 - epsilon)) + A_d_S(1)^((1 - alpha) * (1 - epsilon)))^(1 / ((1 - alpha) * (1 - epsilon))) * (1 - alpha + tau(1) * A_c_S(1)^((1 - alpha) * (1 - epsilon)) / (A_c_S(1)^((1 - alpha) * (1 - epsilon)) + (1 + tau(1))^epsilon * A_d_S(1)^((1 - alpha) * (1 - epsilon))));
% Update environmental quality
S(n) = min(max(0.00000000000000001, -qsi * (Y_d_N(n-1) + Y_d_S(n-1)) + (1 + delta) * S(n-1)), S_bar);
    %Consumption in the base paper is B.6 in the Appendix
%%%Simulations
for n = 2:numsim
A_c_N(n) = (1 + gamma * eta_c_N * s_c_N(n)) * A_c_N(n-1);
A_d_N(n) = (1 + gamma * eta_d_N * (1 - s_c_N(n))) * A_d_N(n-1);
A_c_S(n) = (kk* s_c_S(n) * A_c_N(n)) + ((1- kk)* (1-s_c_S(n))* A_c_S(n-1);
A_d_S(n) = (kk* s_c_S(n) * A_d_N(n)) + ((1- kk)* (1-s_c_S(n))* A_d_S(n-1);
Q(n) = max(0, Q(n-1) - (k2 * alpha^(alpha * (1 - epsilon)) * alpha2 * (cost(Q(n-1)) * (1 + theta(n-1)))^(alpha2 * (1 - epsilon) - 1) * (A_c_N(n)^(1 + phi)) * A_d_N(n)^((1 - alpha1) / (1 - alpha)) / ((k1 * A_d_N(n)^phi1 + ((alpha^alpha) * (1 + tau(n)) * (cost(Q(n-1)) * (1 + theta(n-1)))^alpha2 * A_c_N(n)^(1 - alpha))^(1 - epsilon))^(1 / phi) * (k1 * (1 + tau(n))^epsilon * A_d_N(n)^phi1 + (alpha^alpha * (cost(Q(n-1)) * (1 + theta(n-1)))^alpha2 * A_c_N(n)^(1 - alpha))^(1 - epsilon)))));

% Compute consumption for North and South
%The North will be the same as in the old model, except with some North subscripts. The South should change    
C_N(n) = (alpha / psi)^(alpha / (1 - alpha)) * A_c_N(n) * A_d_N(n) / ((1 + tau(n))^(1 - epsilon) * A_c_N(n)^((1 - alpha) * (1 - epsilon)) + A_d_N(n)^((1 - alpha) * (1 - epsilon)))^(1 / ((1 - alpha) * (1 - epsilon))) * (1 - alpha + tau(n) * A_c_N(n)^((1 - alpha) * (1 - epsilon)) / (A_c_N(n)^((1 - alpha) * (1 - epsilon)) + (1 + tau(n))^epsilon * A_d_N(n)^((1 - alpha) * (1 - epsilon))));
%The South uses other technologies, so the consumption parameter does change. This needs to be accounted for in the South model.
%I have updated the new productivity parameter for the South, so this
%should now be complete
C_S(n) = (alpha / psi)^(alpha / (1 - alpha)) * A_c_S(n) * A_d_S(n) / ((1 + tau(n))^(1 - epsilon) * A_c_S(n)^((1 - alpha) * (1 - epsilon)) + A_d_S(n)^((1 - alpha) * (1 - epsilon)))^(1 / ((1 - alpha) * (1 - epsilon))) * (1 - alpha + tau(n) * A_c_S(n)^((1 - alpha) * (1 - epsilon)) / (A_c_S(n)^((1 - alpha) * (1 - epsilon)) + (1 + tau(n))^epsilon * A_d_S(n)^((1 - alpha) * (1 - epsilon))));

% Update profit ratios for North and South
profitratio_N = eta_c_N * (1 + gamma) * (1 - alpha) * (1 + tau(n))^(1 - alpha) * A_c_N(n-1) / (eta_d_N * (1 + gamma * eta_d_N)^(phi + 1) * A_d_N(n) * (cost(Q(n-1)) * (1 + theta(n-1)))^(alpha2 * (epsilon - 1)) * A_d_N(n)^phi1);
profitratio_S = kappa * (1 - alpha) * (1 + tau(n))^(1 - alpha) * A_c_S(n-1) / (eta_d_S * (1 + gamma * eta_d_S)^(phi + 1) * A_d_S(n) * (cost(Q(n-1)) * (1 + theta(n-1)))^(alpha2 * (epsilon - 1)) * A_d_S(n)^phi1);

% Use these profit ratios for further calculations
% the profit ratios need to be calculated from the new paper. check
% this and add it to the model, look where this has been put in in the old
% the old one

% Update environmental quality
S(n) = min(max(0.00000000000000001, -qsi * (Y_d_N(n-1) + Y_d_S(n-1)) + (1 + delta) * S(n-1)), S_bar);
end

% Compute utility
Teste1 = (1 + repmat(rho, numsim, 1));
Teste2 = Teste1.^((0:numsim-1)');
Teste3 = 1 ./ Teste2;

Util = zeros(1, numsim);
for j = 1:numsim
Util(j) = -(1 / (1 - sigma)) * Teste3(j) * (phiS(S(j)) * C(j))^(1 - sigma);
end

% Compute the clean research subsidy
for k = 1:numsim
if s_c(k) == 0 && profitratio_N(k) > 1
q(k) = 0;
elseif s_c(k) == 1 && profitratio_N(k) < 1
            q(k) = 0;
        else
            q(k) = profitratio_N(k) - 1;
        end
    end

    Resp.Util = Util; % utility flow
    Resp.C = C; % consumption
    Resp.S = S; % quality of the environment
    Resp.Ac = A_c_N; % productivity of the clean sector (North)
    Resp.Ad = A_d_N; % productivity of the dirty sector (North)
    Resp.tau = tau; % input tax
    Resp.S_c = s_c; % share of scientists in clean research
    Resp.Q = Q; % stock of resource
    Resp.theta = theta; % resource tax
    Resp.R = R; % resource extraction
    Resp.Y_d = Y_d_N; % production of dirty input (North)
    Resp.cleansubs = cleansubs; %subsidy to clean research


