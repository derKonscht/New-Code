% This function takes as arguments the paths for the share of scientists in the clean sector and the input tax, and delivers a tax which is the maximum between the tax necessary to implement the allocation of scientists - without any subsidy - and the tax given in the argument.
%the environmental policy does only happen in the North, so the tax is also
%only applied there. %this makes sure that the model's policy implementation is feasible without the need for subsidies
function tax = makeitasNS(s,t)
global numsim phi eta_c_N eta_d_N epsilon gamma Ac0 Ad0
taxx=zeros(numsim,1);
Acdec=Ac0*ones(numsim,1);
Addec=Ad0*ones(numsim,1);
if s(1)==1
       taxx(1) = (eta_d_N/eta_c_N*(1+gamma*eta_d_N)^(-phi-1)*(Acdec(1)/Addec(1))^phi)^(1/epsilon)-1; %this means that if all scientists are allocated to the clean sector, it calculates the necessary tax rate for this
elseif s(1)>0
       taxx(1) = (eta_d_N/eta_c_N*(1+gamma*eta_d_N*(1-s(1)))^(-phi-1)*(1+gamma*eta_c_N*s(1))^(phi+1)*(Acdec(1)/Addec(1))^phi)^(1/epsilon)-1; %if only some are allocated, it calculates the tax rate that accounts for this partial allocation
end
for i=2:numsim
    Acdec(i)=(1+gamma*eta_c_N*s(i-1))*Acdec(i-1); %it updates the productivity levels for clean and dirty sectors based on the previous period's allocation, it then calculates the necessary tax rate for the current allocation, similar to the first period
    Addec(i)=(1+gamma*eta_d_N*(1-s(i-1)))*Addec(i-1);
    if s(i)==1
       taxx(i)=(eta_d_N/eta_c_N*(1+gamma*eta_d_N)^(-phi-1)* (Acdec(i)/Addec(i))^phi)^(1/epsilon)-1;
    elseif s(i)>0
        taxx(i)=(eta_d_N/eta_c_N*(1+gamma*eta_d_N*(1-s(i)))^(-phi-1)*(1+gamma*eta_c*s(i))^(phi+1)*(Acdec(i)/Addec(i))^phi)^(1/epsilon)-1;
    end
end
tax=max(1.0*taxx,t);
end
