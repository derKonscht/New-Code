%Parameter calculation:
%1. new Ac0 and Ad0
global dt rho alpha alpha1 alpha2 eta_d eta_c pf psi phi phi1 kappa k1 k2 qsi epsilon sigma gamma numsim t_disaster
dt = 5;
rho = 0.001*dt; % discount rate
epsilon = 10;  % elasticity of substitution, must be greater than 1 for the program to work properly.
sigma = 2; % Coefficient of relative risk aversion 
alpha = 1/3; % share of machines in production
gamma = 1; % size of innovation
numsim = 80; % number of periods
t_disaster = 6;

emission0 = 23.04;
Yc0 = 318.275;
Yd0 = 2545.639;
R00 = 567.7704867;
Q00 = 5727.801057; % stock of resource in 2018 (in billion of barrels)
Q0 = Q00 - R00; % stock of resource in 2022 (in billion of barrels)
proil0 = 1000000 * (-29.62811308 * Q00 + 0.00466266 * Q00^2 + 80971.56768); % price of resource in 2002 in 2000$
p0c0 = 4.95768E+14; % GDP in 2018 - 2022 in 2017$
spendRgdp = (proil0 * R00) / p0c0; % ratio spending on resource over GDP

alpha2 = (1-alpha) * (1 + (Yd0 / Yc0)^(1/epsilon - 1)) * spendRgdp;
alpha1 = alpha - alpha2;
eta_d = (1-alpha) / (1-alpha1) * 0.02 * dt; % probability of success in the dirty sector
eta_c = 0.02 * dt; % probability of success in the clean sector
psi = alpha^2; % cost of machines (the expressions for Ad0 and Ac0 assume this relationship)
phi = (1-alpha) * (1-epsilon);
phi1 = (1-alpha1) * (1-epsilon);
pf = alpha2^(-1) * proil0 * R00 / Yd0 * (1 + (Yd0 / Yc0)^(1/epsilon - 1))^(1/(1-epsilon));
cos = 1000000 * (-29.62811308 * Q00 + 0.00466266 * Q00^2 + 80971.56768) / pf;
cost_value =1000000*(-394.726* Q00 +0.0421224* Q00 ^2 + 942115)/pf;
kappa = (1-alpha) / (1-alpha1) * (psi^alpha2 * alpha1^alpha1 * alpha2^alpha2 / alpha^alpha)^(1-epsilon);
k1 = (psi^alpha2 * alpha1^alpha1 * alpha2^alpha2)^(1-epsilon);
k2 = (psi^(-alpha1/(1-alpha)) * alpha^(alpha/(1-alpha)) * alpha1^(alpha1/(1-alpha)) * alpha2^(alpha2/(1-alpha)));
qsi = emission0 / Yd0;
for n=2:numsim
pc = (1-((1+tau(n))*pd)^(1-epsilon))^(1/(1-epsilon));
pd = fsolve (@ (pd) psi^(alpha1/(1-alpha1)-alpha/(1-alpha))*alpha^(alpha/(1-alpha))*(1 - pd^(1-epsilon)*(1+tau(n))^(1-epsilon))^(1/phi)*A_c(n)/(alpha1^(alpha1/(1-alpha1))*R(n)^(alpha2/(1-alpha1))*(1+tau(n))^(epsilon*alpha2/(1-alpha1))*(1-tau(n)*(1+tau(n))^(-epsilon)*pd^(1-epsilon))^(alpha2/(1-alpha1))*pd^((1-alpha2*(1-epsilon))/(1-alpha1))*A_d(n)) - 1, pd0);
end
% Calculate initial productivities new
Ac0 = (alpha / psi)^(-alpha / (1-alpha)) * Yc0 * (1 + (Yc0 / Yd0)^(1/epsilon - 1))^(alpha / phi + 1); 
Ad0 = (alpha / psi)^(-alpha / (1-alpha)) * Yd0 * (1 + (Yd0 / Yc0)^(1/epsilon - 1))^(alpha / phi + 1);
Ad0ex = (alpha1^alpha1 * alpha2^alpha2 / (alpha^(2*alpha1)))^(-1 / (1-alpha1)) * Yd0^((1-alpha) / (1-alpha1)) * (1 + (Yd0 / Yc0)^(1/epsilon - 1))^((alpha + phi) / phi1) * cost_value^(alpha2 / (1-alpha1));
Ac0ex = (alpha1^alpha1 * alpha2^alpha2 / alpha^(alpha1 - alpha2) * cost_value^(-alpha2) * Ad0^(1-alpha1) * (Yc0 / Yd0)^(1/epsilon))^(1 / (1-alpha));

% Display calculated values
disp('Ac0: '), disp(Ac0)
disp('Ad0: '), disp(Ad0)
disp('Ad0ex: '), disp(Ad0ex)
disp('Ac0ex: '), disp(Ac0ex)

%This is for the initialising of the simulation loop, now take the 