function U = opttaxNS(x, Ac0, Ad0,S0) 
global rho sigma psi alpha gamma eta_d eta_c qsi epsilon delta numsim S_bar kk
%This is the Utility calculation for the NS case, we initialize the new
%parameters and from there, the loop will run for numsim periods to
%optimize the Utility, which will then be used in the mainprogram.
%   Detailed explanation goes here
%kappa is used in the mainprogramexhaust. So for the profitability of A_c_S
%and the probability of success we introducde a new variable "kk"
%%% Setting vectors' sizes
A_c_N = zeros(numsim,1);
A_d_N = zeros(numsim,1);
A_c_S= zeros(numsim,1);
A_d_S = zeros(numsim,1);
C = zeros(numsim,1);
S = zeros(numsim,1);
s_c = x(1:1:numsim);
tau = x(numsim+1:1:2*numsim);
%%% Initial values %Ac0 and Ad0 are calculated in the mainprogram
%%%Initial values
S(1) = S0; 
A_c_N(1) = (1 + gamma * eta_c_N * s_c_N(1)) * Ac0;
A_d_N(1) = (1 + gamma * eta_d_N * (1 - s_c_N(1))) * Ad0;
A_c_S(1) = (kk* s_c_S(1) * A_c0) + ((1- kk)* (1-s_c_S(1))* Ac0;
A_d_S(1) = (kk* s_c_S(1) * A_c0) + ((1- kk)* (1-s_c_S(1))* Ad0;
%kappa is already used in the mainprogram
%Consumption comes from B.6 in AABH, now do it with the updated
%constraints especially from the South with kappa
C_N(1) = (alpha / psi)^(alpha / (1 - alpha)) * A_c_N(1) * A_d_N(1) / ((1 + tau(1))^(1 - epsilon) * A_c_N(1)^((1 - alpha) * (1 - epsilon)) + A_d_N(1)^((1 - alpha) * (1 - epsilon)))^(1 / ((1 - alpha) * (1 - epsilon))) * (1 - alpha + tau(1) * A_c_N(1)^((1 - alpha) * (1 - epsilon)) / (A_c_N(1)^((1 - alpha) * (1 - epsilon)) + (1 + tau(1))^epsilon * A_d_N(1)^((1 - alpha) * (1 - epsilon))));
C_S(1) = (alpha / psi)^(alpha / (1 - alpha)) * A_c_S(1) * A_d_S(1) / ((1 + tau(1))^(1 - epsilon) * A_c_S(1)^((1 - alpha) * (1 - epsilon)) + A_d_S(1)^((1 - alpha) * (1 - epsilon)))^(1 / ((1 - alpha) * (1 - epsilon))) * (1 - alpha + tau(1) * A_c_S(1)^((1 - alpha) * (1 - epsilon)) / (A_c_S(1)^((1 - alpha) * (1 - epsilon)) + (1 + tau(1))^epsilon * A_d_S(1)^((1 - alpha) * (1 - epsilon))));

%%% Simulations, the North has the same technology as before, the South
%%% cannot innovate on its own, so the average productivity differs and
%%% takes another form
for n = 2:numsim
A_c_N(n)=(1+gamma*eta_c*s_c(n))*A_c_N(n-1);
A_d_N(n)=(1+gamma*eta_d*(1-s_c(n)))*A_d_N(n-1);
A_c_S(n) = (kk* s_c_S(n)* A_c_S(n)+ ((1-kk)* (1- s_c_S(n))* A_c_S(n-1);
A_d_S(n) = (kk* s_c_S(n) * A_d_N(n) + ((1- kk)* (1-s_c_S(n))* A_d_S(n-1);
C_N(n)=(alpha/psi)^(alpha/(1-alpha))*A_c_N(n)*A_d_N(n)/((1+tau(n))^(1-epsilon)*A_c_N(n)^((1-alpha)*(1-epsilon))+A_d_N(n)^((1-alpha)*(1-epsilon)))^(1/((1-alpha)*(1-epsilon)))*(1-alpha+tau(n)*A_c_N(n)^((1-alpha)*(1-epsilon))/(A_c_N(n)^((1-alpha)*(1-epsilon))+(1+tau(n))^epsilon*A_d_N(n)^((1-alpha)*(1-epsilon))));
C_S(n)=(alpha/psi)^(alpha/(1-alpha))*A_c_S(n)*A_d_S(n)/((1+tau(n))^(1-epsilon)*A_c_S(n)^((1-alpha)*(1-epsilon))+A_d_S(n)^((1-alpha)*(1-epsilon)))^(1/((1-alpha)*(1-epsilon)))*(1-alpha+tau(n)*A_c_S(n)^((1-alpha)*(1-epsilon))/(A_c_S(n)^((1-alpha)*(1-epsilon))+(1+tau(n))^epsilon*A_d_S(n)^((1-alpha)*(1-epsilon))));

%then this is the new environemtnal constraint in the North South model
%Y_d and Y_c for the North and the South
  S(n) = min(max(0.00000000000000001, -qsi * (Y_d_N(n-1) + Y_d_S(n-1)) + (1 + delta) * S(n-1)), S_bar);% S should not be very close to zero
end
% Assuming C_N and C_S are vectors containing consumption values for North and South
% And assuming phiS(S) is the environmental quality factor for the combined environment

% Discount factor initialization
Teste1 = (1 + repmat(rho, numsim, 1));
Teste2 = Teste1.^((0:numsim-1)');
Teste3 = 1 ./ Teste2;

% Utility calculation
Util = zeros(numsim, 1);

for j = 1:numsim
    % Combined utility for North and South
    Util(j) = -(1 / (1 - sigma)) * Teste3(j) * (phiS(S(j)) * (C_N(j)*C_S(j)))^(1 - sigma);
end
% Total utility (sum of discounted utilities)
U = sum(Util);



