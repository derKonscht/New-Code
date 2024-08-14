%here all the changes in the code are written down, to get an overview on
%what we need to adjust and what can stay the same

%1. Environmental constraint
S(n) = min(max(0.00000000000000001, -qsi * (Y_d_N(n-1) + Y_d_S(n-1)) + (1 + delta) * S(n-1)), S_bar);

%2. profit ratios change

profit_N = eta_c_N * (1 + gamma) * (1 - alpha) * (1 + tau(n))^(1 - alpha) * A_c_N(n-1) / (eta_d_N * (1 + gamma * eta_d_N)^(phi + 1) * A_d_N(n) * (cost(Q(n-1)) * (1 + theta(n-1)))^(alpha2 * (epsilon - 1)) * A_d_N(n)^phi1);
profit_S = kappa_j * (1 - alpha) * (1 + tau(n))^(1 - alpha) * A_c_S(n-1) / (eta_d_S * (1 + gamma * eta_d_S)^(phi + 1) * A_d_S(n) * (cost(Q(n-1)) * (1 + theta(n-1)))^(alpha2 * (epsilon - 1)) * A_d_S(n)^phi1);

%3. relative profit from undertaking research in a sector
 % Update profit ratios for North and South
 profitratio_N = (eta_c_N / eta_d_N) * (1 / (1 - alpha)) * (A_c_N(n-1) / A_d_N(n-1));
 profitratio_S = (kappa_c / kappa_d) * (1 / (1 - alpha)) * (A_c_S(n) / A_N(n))^(-phi) * (A_N(n) / A_d_N(n));
 
 %4. Define the new profitability A and C for both North and South
 A_c_N(1) = (1 + gamma * eta_c_N * s_c(1)) * Ac0; %this is the same as in the original paper: equation 11 
 A_d_N(1) = (1 + gamma * eta_d_N * (1 - s_c(1))) * Ad0; %equation 11 can be found in the mysimenvtaxnew2 code part
 A_c_S(n) = kappa_c * s_c(n) * A_c_N(n) + (1 - kappa_c * s_c(n)) * A_c_S(n-1); %this needs to be adjusted to the one in the North South model
 A_d_S(n) = kappa_d * (1 - s_c(n)) * A_d_N(n) + (1 - kappa_d * (1 - s_c(n))) * A_d_S(n-1); %now we have described the evolution of productivity in the South and in the North
 C_N(1) = (alpha / psi)^(alpha / (1 - alpha)) * A_c_N(1) * A_d_N(1) / ((1 + tau(1))^(1 - epsilon) * A_c_N(1)^((1 - alpha) * (1 - epsilon)) + A_d_N(1)^((1 - alpha) * (1 - epsilon)))^(1 / ((1 - alpha) * (1 - epsilon))) * (1 - alpha + tau(1) * A_c_N(1)^((1 - alpha) * (1 - epsilon)) / (A_c_N(1)^((1 - alpha) * (1 - epsilon)) + (1 + tau(1))^epsilon * A_d_N(1)^((1 - alpha) * (1 - epsilon))));
 C_S(1) = (alpha / psi)^(alpha / (1 - alpha)) * A_c_S(1) * A_d_S(1) / ((1 + tau(1))^(1 - epsilon) * A_c_S(1)^((1 - alpha) * (1 - epsilon)) + A_d_S(1)^((1 - alpha) * (1 - epsilon)))^(1 / ((1 - alpha) * (1 - epsilon))) * (1 - alpha + tau(1) * A_c_S(1)^((1 - alpha) * (1 - epsilon)) / (A_c_S(1)^((1 - alpha) * (1 - epsilon)) + (1 + tau(1))^epsilon * A_d_S(1)^((1 - alpha) * (1 - epsilon))));
 %There is another C if we take the exhaustible variable case into account
 C(1) = k2*A_c(1)*A_d(1)^((1-alpha1)/(1-alpha))*((1-alpha)*(1+tau(1))^epsilon*k1*A_d(1)^phi1+(1+tau(1)-(alpha1+alpha2/(1+theta(1))))*(alpha^alpha*(cost(Q(1))*(1+theta(1)))^alpha2*A_c(1)^(1-alpha))^(1-epsilon))/(D1*D2);
 %5. the important thing is the to check that the initial values are
 %calculated differently. e.g the Starting point of A_c0 and A_d0 as shown
 %in the mainprogramexhaustible case might change in the North South model.
 %Try to calculate the new starting point 
 
 %6. try different utility functions by changing mysimenvtaxnew2, maybe the
 %utility function found in the Early Version of the paper
 % Compute utility using the new function
    Util = zeros(1, numsim);
    for j = 1:numsim
        min_S = min(S(j), S_bar);
        term = 2 * sqrt(5 * min_S) - min_S;
        Util(j) = ((term * (C_N(j) + C_S(j)))^(1 - 1.4)) / (1 - 1.4);
    end

    % Compute discounted utility
    Teste1 = (1 + repmat(rho, numsim, 1));
    Teste2 = Teste1.^((0:numsim-1)');
    Teste3 = 1 ./ Teste2;
    discountedUtil = Teste3' .* Util;

    % Summing the discounted utility
    totalUtil = sum(discountedUtil);

%7. new variables and parameter calculations   
%some parameters are calculated in the calibrationexhaustible script. This
%will be the parameters for the North calculated that way
pf=alpha2^(-1)*proil0*R00/Yd0*(1+(Yd0/Yc0)^(1/epsilon-1))^(1/(1-epsilon));
Ad0 = (alpha1^alpha1*alpha2^alpha2/(alpha^(2*alpha1)))^(-1/(1-alpha1))*Yd0^((1-alpha)/(1-alpha1))*(1+(Yd0/Yc0)^(1/epsilon-1))^((alpha+phi)/phi1)*cost(Q00)^(alpha2/(1-alpha1));
Ac0=(alpha1^alpha1*alpha2^alpha2/alpha^(alpha1-alpha2)*cost(Q00)^(-alpha2)*Ad0^(1-alpha1)*(Yc0/Yd0)^(1/epsilon))^(1/(1-alpha));
alpha1= alpha- alpha2;
alpha2= (1-alpha)*(1+(Yd0/Yc0)^((1/epsilon)-1))*spendRgdp;

% 8. Parameters linked to the resource: The numbers come from the new
% calculations done for this. 
R00=567.7704867; % resource production in 2018 - 2022 (in billions of barrels)
Q00=5727.801057; % stock of resource in 2018 (in billion of barrels)
Q0=Q00-R00; % stock of resource in 2022 (in billion of barrels)
proil0 = 1000000*(-29.62811308*Q00 +0.00466266*Q00^2 + 80971.56768); % price of resource in 2002 in 2000$
p0c0= 4.95768E+14; % GDP in 2018 - 2022 in 2017$
spendRgdp= (proil0*R00)/p0c0; % ratio spending on resource over GDP
eta_d = (1-alpha)/(1-alpha1)*0.02*dt; % probability of success in the dirty sector
eta_c = 0.02*dt; % probability of success in the the clean sector (Check whether these are different in the North and South