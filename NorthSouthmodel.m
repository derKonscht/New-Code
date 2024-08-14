%This is the calculation part for the North South model. I introduce new
%parameters for the South and the North respectively. After that running
%the mainprogram yields then optimal tax path for the North and the South
global dt numsim alpha_N alpha_S rho epsilon qsi delta gamma eta_cN eta_dN eta_cS eta_dS kappa_c kappa_d psi phi phi1

% Parameters
dt = 5;
numsim = 80; %time periods
alpha_N = 0.33; % Example parameter for North
alpha_S = 0.33; % Example parameter for South
rho = 0.01*dt; % Discount rate (either=0.001 or 0.015)
epsilon = 1.5; % Elasticity of substitution
delta = 0.01; % Rate of environmental regeneration
gamma = 0.1; % Innovation step size
eta_cN = 0.8; % Innovation probability for clean tech in North
eta_dN = 0.8; % Innovation probability for dirty tech in North
eta_cS = 0.5; % Imitation probability for clean tech in South
eta_dS = 0.5; % Imitation probability for dirty tech in South
kappa_c = 0.6; % Imitation factor for clean tech in South
kappa_d = 0.6; % Imitation factor for dirty tech in South
qsi = emission0/Yd0; %rate of environmental degradation
emission0= 23.04; %current emissions in ppm
% Initial conditions
S_bar = 1; % Maximum possible environmental quality
S_N0 = S_bar; % Initial environmental quality for North
S_S0 = S_bar; % Initial environmental quality for South


% Initialize variables
S_N = zeros(1, numsim);
S_S = zeros(1, numsim);
S_N(1) = S_N0;
S_S(1) = S_S0;
Y_N = zeros(1, numsim);
Y_S = zeros(1, numsim);
C_N = zeros(1, numsim);
C_S = zeros(1, numsim);
A_cN = zeros(1, numsim);
A_dN = zeros(1, numsim);
A_cS = zeros(1, numsim);
A_dS = zeros(1, numsim);
% 5. preliminary computations
psi= alpha^2; % cost of machines (the expressions for Ad0 and Ac0 assume this relationship)
phi=(1-alpha)*(1-epsilon);
phi1=(1-alpha1)*(1-epsilon);
kappa = (1-alpha)/(1-alpha1)*(psi^alpha2*alpha1^alpha1*alpha2^alpha2/alpha^alpha)^(1-epsilon);
k1 = ((psi^alpha2)*(alpha1^alpha1)*(alpha2^alpha2))^(1-epsilon);
k2 = (psi^(-alpha1/(1-alpha)))*(alpha^(alpha/(1-alpha)))*(alpha1^(alpha1/(1-alpha)))*(alpha2^(alpha2/(1-alpha)));
pf=alpha2^(-1)*proil0*R00/Yd0*(1+(Yd0/Yc0)^(1/epsilon-1))^(1/(1-epsilon));
Ad0 = (alpha1^alpha1*alpha2^alpha2/(alpha^(2*alpha1)))^(-1/(1-alpha1))*Yd0^((1-alpha)/(1-alpha1))*(1+(Yd0/Yc0)^(1/epsilon-1))^((alpha+phi)/phi1)*cost(Q00)^(alpha2/(1-alpha1));
Ac0=(alpha1^alpha1*alpha2^alpha2/alpha^(alpha1-alpha2)*cost(Q00)^(-alpha2)*Ad0^(1-alpha1)*(Yc0/Yd0)^(1/epsilon))^(1/(1-alpha));
qsi = emission0/Yd0;
D1= (k1*A_d(1)^phi1+((alpha^alpha)*(1+tau(1))*(cost(Q(1))*(1+theta(1)))^alpha2*A_c(1)^(1-alpha))^(1-epsilon))^(1/phi);
D1 = (k1*A_d(n)^phi1+((alpha^alpha)*(1+tau(n))*(cost(Q(n))*(1+theta(n)))^alpha2*A_c(n)^(1-alpha))^(1-epsilon))^(1/phi);

%********************************************************************
%********************************************************************
%2. options
mode = 1; %1: for the model with directed technical change and both instruments, 0.5: for a model with DTC but only the carbon tax, 0: for a model without dtc
% the program for the case with a carbon tax only is correct only if epsilon>(2-alpha)/(1-alpha), moreover if this is the case, it assumes that when several allocation of scientists are possible equilibria, the equilibrium chosen is the interior one.
delay = 0; % delay before applying optimal policy (possible only for mode = 1)
display_iter=0;    %=1 if you want to see iterations, ~1 if otherwise
display_diag=1;    %=1 if you want to see diagnostics, ~1 if otherwise


%********************************************************************
%********************************************************************
%3. initial guesses
if mode == 1
    sc0=sce10rho01; % initial guess for the share of scientists in clean research
    t0=taue10rho01; % initial guess for the input tax
    x0=[sc0; t0]; % vector stacking the guess for the share of scientists in clean research and the guess for the input tax
    ub=[ones(numsim,1); 100000*ones(numsim,1)]; % upper bound for the optimization
    lb=[zeros(numsim,1); zeros(numsim,1)]; % lower boud for the optimization
else
    t0=taue10rho01taxonly; % initial guess for the input tax
    x0=t0;
    ub=100000*ones(numsim,1);
    lb=zeros(numsim,1);
end

%**** No intervention in the code is necessary beyond this point

%********************************************************************
%********************************************************************
% 4. Warnings
if epsilon <= 1
    display ('*********the program is not stable for that range******************')
end
if  mode ==.5 && epsilon <= (2-alpha)/(1-alpha)
    display('*********epsilon is too low for the code to make sense ******************')
end
if  mode <1 && delay>0
    display('*********delay is a possible option only with DTC and two instruments******************')
end


%********************************************************************
%********************************************************************
% 5. preliminary computations
phi=(1-alpha)*(1-epsilon);
qsi = emission0/Yd0;
Ac0=(alpha/psi)^(-alpha/(1-alpha))*Yc0*(1+(Yc0/Yd0)^(1/epsilon-1))^(alpha/phi+1); % corresponding initial value for Ac0 (assuming that the optimal subsidy to the use of all machines is implemented)
Ad0 =(alpha/psi)^(-alpha/(1-alpha))*Yd0*(1+(Yd0/Yc0)^(1/epsilon-1))^(alpha/phi+1); % and for Ad0
numsimleft=numsim-delay;
StockAc=[];
StockAd=[];
StockC=[];
StockQ=[];
Stocktau=[];
Stocksc=[];
StockS=[];
utcomp=0;
% evolution of the economy during the delay period
numsim=numsimleft;

%********************************************************************
%********************************************************************
%6. optimization
if display_iter==1
    if display_diag==1
    options = optimset('Display','iter', 'Diagnostics','on','FunValCheck','on','TolFun',1e-11,'TolX',1e-9);
    else
    options = optimset('Display','iter', 'Diagnostics','off', 'FunValCheck','on','TolFun',1e-11,'TolX',1e-9);
    end
else
    if display_diag==1
    options = optimset('Display','off', 'Diagnostics','on','FunValCheck','on','TolFun',1e-11,'TolX',1e-9);
    else 
    options = optimset('Display','off', 'Diagnostics','off', 'FunValCheck','on','TolFun',1e-11,'TolX',1e-9);
    end
end    
if mode == 1
    [x,fval,exitflag] = fmincon(@(x)mysimopttaxnew2NorthSouth(x, Ac0, Ad0, S0),x0,[],[],[],[],lb,ub,[],options); % optimization with DTC and 2 instruments
    % mysimopttaxnew2 computes the utility given x - a vector stacking the path for the share of scientists in clean research and the input tax - and the initial values of clean and dirty productivities and quality of the environment
    else
    if mode == 0
        [x,fval,exitflag] = fmincon(@(x)mysimopttaxnew2noDTCNorthSouth(x, Ac0, Ad0, S0),x0,[],[],[],[],lb,ub,[],options); % optimization without DTC
        % similarly mysimopttaxnew2noDTC computes the utility given a vector of input tax and assuming that there is no DTC
        x=[eta_d/(eta_c+eta_d)*ones(numsim,1);x]; % x now combines the (constant) allocation of scientists with the tax rate
    else
        [x,fval,exitflag] = fmincon(@(x)mysimopttaxnew2onlytauNorthSouth(x, Ac0, Ad0, S0),x0,[],[],[],[],lb,ub,[],options); % optimization with DTC and carbon tax only
        % mysimopttaxnew2onlytau computes the utility given a vector of input tax x and deriving the corresponding allocation for scientists
        % if several allocation of scientists are equilibria, the interior one is chosen
    end
end
if exitflag<=0
    display('*********problem, optimization not converged******************')
end


%********************************************************************
%********************************************************************
%7. take results
Util=utcomp-1/(1+rho)^delay*fval;
Resp = mysimenvtaxnew2NorthSouth(x, Ac0, Ad0, S0);
    % mysimenvtaxnew2 computes all the relevant parameters of the economy given the vector x which stacks the share of scientists in clean research and the input tax
    xx=[Stocksc;x(1:numsim)];
    Acc=[StockAc;Resp.Ac];
    Add=[StockAd;Resp.Ad];
    Qq=[StockQ;Resp.Q]; % Qq is the clean research subsidy, it makes sense only if mode=1
    tauu=[Stocktau;Resp.tau];
    Cc=[StockC;Resp.C];
    Ss=[StockS;Resp.S];
Tt=zeros(numsim+delay,1);
Ratio=zeros(numsim+delay,1);
for i=1:(numsim+delay)
    Tt(i)=3/log(2)*log((280*2^(t_disaster/3)-Ss(i))/280);
    Ratio(i)=1/(1+(Add(i)/Acc(i))^(epsilon*(1-alpha))*(1+tauu(i))^(-epsilon));
end


%********************************************************************
%********************************************************************
%8. plot results

numsimc=numsim+delay;
subplot(2,2,1)
plot(linspace(1,numsimc,numsimc),[xx(1:numsimc),ones(numsimc,1)-xx(1:numsimc)]);
xlabel('period')
ylabel('Fraction of Scientists')
legend('S_c','S_d')
subplot(2,2,2)
plot(linspace(1,numsimc,numsimc),[Acc(1:numsimc),Add(1:numsimc),Cc(1:numsimc)]);
xlabel('period')
ylabel('economy')
legend('Ac','Ad','C')
subplot(2,2,3)
plot(linspace(1,numsimc,numsimc),tauu(1:numsimc));
xlabel('period')
ylabel('input tax')
legend('tau')
subplot(2,2,4)
plot(linspace(1,numsimc,numsimc),Tt(1:numsimc));
xlabel('period')
ylabel('temp')
legend('delta')
