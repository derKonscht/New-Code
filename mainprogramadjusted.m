%1. initial parameters
global dt rho sigma psi alpha alpha1 gamma eta_d eta_c qsi epsilon delta numsim phi S_bar lambda t_disaster max_t 
dt = 5;  % number of years in a period
rho = 0.001 * dt; % discount rate
epsilon = 10;  % elasticity of substitution, must be greater than 1 for the program to work properly.
sigma = 2; % Coefficient of relative risk aversion 
alpha = 0.33; % share of machines in production
psi = alpha^2; % cost of machines
gamma = 1; % size of innovation
eta_d = [0.02 * dt, (1-alpha)/(1-alpha1)*0.02 * dt]; % probability of success in the dirty sector for North and South
eta_c = [0.02 * dt, 0.02 * dt]; % probability of success in the clean sector for North and South
numsim = 80; % number of periods
qsi= emission0/Yd0; %emission factor from dirty goods
phi= (1-alpha)*(1-epsilon); %per definition

%Now we define the parameters that are North specific and the ones that
%are South specific:
%1. North specific
eta_d_N = 0.02*dt; %probability of success dirty in the North
eta_c_N = 0.02*dt; %probability of success clean in the North
alpha_N= 0.33; 
A_c_N(1) = (1 + gamma * eta_c_N * s_c(1)) * Ac0; %this is the same as in the original paper: equation 11 
A_d_N(1) = (1 + gamma * eta_d_N * (1 - s_c(1))) * Ad0; %equation 11 can be found in the mysimenvtaxnew2 code part
profit_N = eta_c_N * (1 + gamma) * (1 - alpha) * (1 + tau(n))^(1 - alpha) * A_c_N(n-1) / (eta_d_N * (1 + gamma * eta_d_N)^(phi + 1) * A_d_N(n) * (cost(Q(n-1)) * (1 + theta(n-1)))^(alpha2 * (epsilon - 1)) * A_d_N(n)^phi1);
%2. South specific
eta_d_S = 0.02*dt; %probability of success dirty in the South
eta_c_S = 0.02*dt; %probability of success clean in the South
alpha_S = 0.33;
A_c_S(n) = kappa_c * s_c(n) * A_c_N(n) + (1 - kappa_c * s_c(n)) * A_c_S(n-1); %this needs to be adjusted to the one in the North South model
A_d_S(n) = kappa_d * (1 - s_c(n)) * A_d_N(n) + (1 - kappa_d * (1 - s_c(n))) * A_d_S(n-1); %now we have described the evolution of productivity in the South and in the North
profit_S = kappa_j * (1 - alpha) * (1 + tau(n))^(1 - alpha) * A_c_S(n-1) / (eta_d_S * (1 + gamma * eta_d_S)^(phi + 1) * A_d_S(n) * (cost(Q(n-1)) * (1 + theta(n-1)))^(alpha2 * (epsilon - 1)) * A_d_S(n)^phi1);


% initial values for the clean and dirty inputs, source table 11.1 in the Annual Energy Review 2008
Yc0 = [318.26, 318.26]; % production of non-fossil fuel energy in the world primary supply of energy from 2018 to 2022 source in Quadrillion of Btu 
Yd0 = [2545.639, 2545.639]; % production of fossil fuel energy in the world primary supply of energy from 2018 to 2022 source in Quadrillion of Btu
% initial emissions, IAE Data:https://www.iea.org/reports/co2-emissions-in-2022 
emission0 = [17.48251782, 17.48251782]; % world emissions in ppm from 2018 to 2022 (converted from original data IEA in GtCO2)
t_disaster = 6; % disaster temperature
S_bar = 280 * (2^(t_disaster / 3) - 1); % corresponding value for S_bar
S0 = S_bar - 99; % 99ppm increase in CO2 concentration since preindustrial times
max_t = 3; % temperature up to which the damage function is matched with Nordhaus's one
lambda = fmincon(@(l)damage_calibS(l), 0.35, [], [], [], [], 0.00001, 0.999999, [], optimset('Tolfun', 1e-11));
% the formula for Nordhaus's model is used in the function damage_calibS and can be found at http://nordhaus.econ.yale.edu/RICEmodels.htm
delta = 0.5 * emission0 ./ S0; % half of CO2 emissions are absorbed at current atmospheric levels
% Parameters linked to the resource: see Appendix A
R00=567.7704867; % resource production in 2018 - 2022 (in billions of barrels)
Q00=5727.801057; % stock of resource in 2018 (in billion of barrels)
Q0=Q00-R00; % stock of resource in 2022 (in billion of barrels)
proil0 = 1000000*(-29.62811308*Q00 +0.00466266*Q00^2 + 80971.56768); % price of resource in 2022 in 2017$
p0c0= 4.95768E+14; % GDP in 2018 - 2022 in 2017$
spendRgdp= (proil0*R00)/p0c0; % ratio spending on resource over GDP

% Initialize variables
S = zeros(2, numsim);
Q = zeros(2, numsim);
Y = zeros(2, numsim);
C = zeros(2, numsim);
A_c = zeros(2, numsim);
A_d = zeros(2, numsim);
tauu = zeros(2, numsim);
Ratio = zeros(2, numsim);

% Set initial values
S(:, 1) = S0;

% Production function for both countries
for t = 1:numsim
    for k = 1:2 % Loop over North (k=1) and South (k=2)
        Y_c(k, t) = alpha * integral(@(i) A_c(k, t) .* x_c(k, t).^(1 - alpha), 0, 1);
        Y_d(k, t) = (1 - alpha) * integral(@(i) A_d(k, t) .* x_d(k, t).^(1 - alpha), 0, 1);
        Y(k, t) = (Y_c(k, t)^((epsilon - 1) / epsilon) + Y_d(k, t)^((epsilon - 1) / epsilon))^(epsilon / (epsilon - 1));
    end
end

% Environmental quality dynamics
for t = 1:numsim - 1
    for k = 1:2 % Loop over North (k=1) and South (k=2)
        S(k, t + 1) = min(max(S(k, t) - qsi * Y_d(k, t) + delta(k) * S(k, t), 0), S_bar);
    end
end
% Profit maximization
for t = 1:numsim
    for k = 1:2 % Loop over North (k=1) and South (k=2)
        % Define the Lagrangian
        Lagrangian = @(x) p(k, t) * L(k, t)^(1 - alpha) * integral(@(i) A(k, t) * x(i)^(1 - alpha), 0, 1) - w(k, t) * L(k, t) - integral(@(i) p(k, t) * x(i), 0, 1);
        % Solve first order conditions
        options = optimoptions('fsolve', 'Display', 'none');
        x_opt = fsolve(Lagrangian, x0, options);
        x(k, t) = x_opt;
    end
end
for t = 1:numsim - 1
    for k = 1:2 % Loop over North (k=1) and South (k=2)
        % Update production functions
        Y_c(k, t) = alpha * integral(@(i) A_c(k, t) .* x_c(k, t).^(1 - alpha), 0, 1);
        Y_d(k, t) = (1 - alpha) * integral(@(i) A_d(k, t) .* x_d(k, t).^(1 - alpha), 0, 1);
        Y(k, t) = (Y_c(k, t)^((epsilon - 1) / epsilon) + Y_d(k, t)^((epsilon - 1) / epsilon))^(epsilon / (epsilon - 1));
        
        % Update environmental quality
        S(k, t + 1) = min(max(S(k, t) - xi * Y_d(k, t) + delta(k) * S(k, t), 0), S_bar);
        
        % Solve profit maximization and update x
        Lagrangian = @(x) p(k, t) * L(k, t)^(1 - alpha) * integral(@(i) A(k, t) * x(i)^(1 - alpha), 0, 1) - w(k, t) * L(k, t) - integral(@(i) p(k, t) * x(i), 0, 1);
        x_opt = fsolve(Lagrangian, x0, options);
        x(k, t) = x_opt;
        
        % Update consumption
        C(k, t) = Y(k, t) - psi * (integral(@(i) x_c(k, t), 0, 1) + integral(@(i) x_d(k, t), 0, 1)) - c(Q(k, t)) * R(k, t);
        
        % Update innovation
        A_c(k, t + 1) = (1 + eta_c(k) * gamma) * A_c(k, t);
        A_d(k, t + 1) = (1 + eta_d(k) * gamma) * A_d(k, t);
    end
end
%Change the model for the optimization here
if mode == 1
    [x, fval, exitflag] = fmincon(@(x) mysimopttaxnew2(x, Ac0, Ad0, S0), x0, [], [], [], [], lb, ub, [], options); % optimization with DTC and 2 instruments
else
    if mode == 0
        [x, fval, exitflag] = fmincon(@(x) mysimopttaxnew2noDTC(x, Ac0, Ad0, S0), x0, [], [], [], [], lb, ub, [], options); % optimization without DTC
        x = [eta_d(1)/(eta_c(1) + eta_d(1)) * ones(numsim, 1); x]; % x now combines the (constant) allocation of scientists with the tax rate
    else
        [x, fval, exitflag] = fmincon(@(x) mysimopttaxnew2onlytau(x, Ac0, Ad0, S0), x0, [], [], [];

