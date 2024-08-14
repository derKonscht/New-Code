% This is the f function of lemma 1 with an input tax and subsidy differences for the North and South regions.
% s_N and s_S are the shares of scientists in the clean sector in the North and South respectively,
% t_n and t_s are the input taxes, ac_n, ad_n are the productivity levels in the previous period for the North,
% and ac_s, ad_s are the productivity levels in the previous period for the South.
function f_ns = fintNS(s_n, s_s, t_n, t_s, ac_n, ad_n, ac_s, ad_s)
    global eta_c eta_d phi epsilon gamma q_N q_S
    
    % Clean-to-dirty profit ratio in the North
    f_N = (1 + q_N) * (eta_c / eta_d) * (1 + t_n)^epsilon * ...
        ((1 + gamma * eta_c * s_n) / (1 + gamma * eta_d * (1 - s_n)))^(-phi - 1) * (ac_n / ad_n)^(-phi);
    
    % Clean-to-dirty profit ratio in the South, accounting for spillovers
    f_S = (1 + q_S) * (eta_c / eta_d) * (1 + t_s)^epsilon * ...
        ((1 + gamma * eta_c * s_s) / (1 + gamma * eta_d * (1 - s_s)))^(-phi - 1) * (ac_s / ad_s)^(-phi);
    
    % Combine both regions' functions
    f_ns = f_N + f_S; 
end


