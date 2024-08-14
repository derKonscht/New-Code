function U = calculations(S0,AdN0,AcN0,AcS0, AdS0,x)
%This calculates all equations that are adjusted in the North South model
global numsim
%1. Vector sizes:
A_cN = zeros(numsim,1);
A_dN = zeros(numsim,1);
A_cS = zeros(numsim,1);
A_dS = zeros (numsim,1);
C = zeros(numsim,1);
S = zeros(numsim,1);
s_c = x(1:1:numsim);
tau = x(numsim+1:1:2*numsim);

%2. Setting initial values
S(1)= S0; %initial environment constraint
A_cN(1)=(1+gamma*eta_c*s_c(1))*AcN0; %clean production in the North
A_dN(1)=(1+gamma*eta_d*(1-s_c(1)))*AdN0;%dirty prodcution in the South
A_cS(1)= kappa*s_cS(1)*AcN0+(1-kappa*(1-s_cS(1))*AcS0; %initial clean imitation in the South
A_dS(1)= kappa*s_dS(1)*AdN0+(1-kappa*(1-s_dS(1))*AdS0; %initial dirty input production in the South (maybe this is wrong as why would the South would want to imitate dirty technologies
C(1)=(alpha/psi)^(alpha/(1-alpha))*A_c(1)*A_d(1)/((1+tau(1))^(1-epsilon)*A_c(1)^((1-alpha)*(1-epsilon))+A_d(1)^((1-alpha)*(1-epsilon)))^(1/((1-alpha)*(1-epsilon)))*(1-alpha+tau(1)*A_c(1)^((1-alpha)*(1-epsilon))/(A_c(1)^((1-alpha)*(1-epsilon))+(1+tau(1))^epsilon*A_d(1)^((1-alpha)*(1-epsilon))));
%3. simulation for n= 2:numsim periods
for n= 2:numsim
A_cN(n)=(1+gamma*eta_c*s_c(n))*AcN(n-1); 
A_dN(n)=(1+gamma*eta_d*(1-s_c(n)))*AdN(n-1);
A_cS(n)= kappa*s_cS(n)*AcN(n-1)+(1-kappa*(1-s_cS(n))*AcS(n-1);
A_dS(n)= kappa*s_dS(n)*AdN(n-1)+(1-kappa*(1-s_dS(n))*AdS(n-1);  
C(1)=(alpha/psi)^(alpha/(1-alpha))*A_c(n)*A_d(n)/((1+tau(n))^(1-epsilon)*A_c(n)^((1-alpha)*(1-epsilon))+A_d(n)^((1-alpha)*(1-epsilon)))^(1/((1-alpha)*(1-epsilon)))*(1-alpha+tau(n)*A_c(n)^((1-alpha)*(1-epsilon))/(A_c(n)^((1-alpha)*(1-epsilon))+(1+tau(n))^epsilon*A_d(n)^((1-alpha)*(1-epsilon))));
end
%now the Utility calculations
Teste1 = (1+repmat(rho,numsim,1));
Teste2 = Teste1.^((0:numsim-1)');
Teste3 = 1./Teste2;
 
Util=zeros(numsim,1);
for j=1:numsim
   Util(j)=-(1/(1-sigma))*Teste3(j)*(phiS(S(j))*C(j))^(1-sigma);
end
  
U=sum(Util(1:end));

