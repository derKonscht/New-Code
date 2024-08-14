% This function computes the utility for the base case depending on x, a vector stacking the paths for the allocation of scientists and the input tax, the initial values for the productivity of the clean and dirty sector and the inital environmental quality.
function U = mysimopttaxnew2(x, Ac0, Ad0, S0)
global rho sigma psi alpha gamma eta_d eta_c qsi epsilon delta numsim S_bar
%%% Setting vectors' sizes
A_c = zeros(numsim,1);
A_d = zeros(numsim,1);
C = zeros(numsim,1);
S = zeros(numsim,1);
s_c = x(1:1:numsim);
tau = x(numsim+1:1:2*numsim);


%%% Initial values
S(1) = S0;          
A_c(1)=(1+gamma*eta_c*s_c(1))*Ac0;
A_d(1)=(1+gamma*eta_d*(1-s_c(1)))*Ad0;
C(1)=(alpha/psi)^(alpha/(1-alpha))*A_c(1)*A_d(1)/((1+tau(1))^(1-epsilon)*A_c(1)^((1-alpha)*(1-epsilon))+A_d(1)^((1-alpha)*(1-epsilon)))^(1/((1-alpha)*(1-epsilon)))*(1-alpha+tau(1)*A_c(1)^((1-alpha)*(1-epsilon))/(A_c(1)^((1-alpha)*(1-epsilon))+(1+tau(1))^epsilon*A_d(1)^((1-alpha)*(1-epsilon))));
%this is from Appendix B.6 from The Environment and DTC
%%% Simulations
for n = 2:numsim
A_c(n)=(1+gamma*eta_c*s_c(n))*A_c(n-1);
A_d(n)=(1+gamma*eta_d*(1-s_c(n)))*A_d(n-1);
C(n)=(alpha/psi)^(alpha/(1-alpha))*A_c(n)*A_d(n)/((1+tau(n))^(1-epsilon)*A_c(n)^((1-alpha)*(1-epsilon))+A_d(n)^((1-alpha)*(1-epsilon)))^(1/((1-alpha)*(1-epsilon)))*(1-alpha+tau(n)*A_c(n)^((1-alpha)*(1-epsilon))/(A_c(n)^((1-alpha)*(1-epsilon))+(1+tau(n))^epsilon*A_d(n)^((1-alpha)*(1-epsilon))));
%also from the environment and DTC Appendix B.6
  S(n) = min(max(0.00000000000000001,-qsi*(alpha/psi)^(alpha/(1-alpha))*(A_c(n-1)^(((1-alpha)*(1-epsilon))+alpha)*A_d(n-1))/(((1+tau(n-1))^(1-epsilon)*A_c(n-1)^((1-alpha)*(1-epsilon))+A_d(n-1)^((1-alpha)*(1-epsilon)))^(alpha/((1-alpha)*(1-epsilon)))*(A_c(n-1)^((1-alpha)*(1-epsilon))+(1+tau(n-1))^epsilon*A_d(n-1)^((1-alpha)*(1-epsilon)))) + (1+delta)*S(n-1)),S_bar); %notice that S should not be very close to zero
%Debugging statements
%disp(['Iteration ', num2str(n)]);
%disp(['A_c(', num2str(n), '): ', num2str(A_c(n))]);
%disp(['A_d(', num2str(n), '): ', num2str(A_d(n))]);
%disp(['C(', num2str(n), '): ', num2str(C(n))]);
end
Teste1 = (1+repmat(rho,numsim,1));
Teste2 = Teste1.^((0:numsim-1)');
Teste3 = 1./Teste2;
 
Util=zeros(numsim,1);

for j=1:numsim
%disp(['Iterationj: ', num2str(j)]);
%disp(['S(n): ', num2str(S(n))]);
%disp(['C(n): ', num2str(C(n))]);
%disp(['Teste3(j): ', num2str(Teste3(j))]); 
   current_phiS = phiS(S(j));
   current_C = C(j);
Util(j) = -(1/(1-sigma)) * Teste3(j) * (current_phiS * current_C)^(1-sigma);
   % Debugging: Print values to trace issues. The utility turns INF at
   % iteration 56. WHY??
   %fprintf('Iteration %d: S = %f, phiS(S) = %f, C = %f, Util = %f\n', j, S(j), current_phiS, current_C, Util(j));
end
U = sum(Util(1:end));
    %disp(['Final Utility: ', num2str(U)]);




