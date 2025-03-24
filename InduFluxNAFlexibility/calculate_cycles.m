%% Code for calculating the equivalent cycles of operation
% Presudo code is described in Paper titled: Long-Term Revenue Estimation for Battery Performing Arbitrage and Ancillary Services
% input the SoC (belong to 0 and 1) trajectory x_adj
% input profit,

function [total_cycle_100] = calculate_cycles(x_adj)

a_minus =0;
a_plus =0;

vec=[];

i=1;
N=length(x_adj);

while i< length(x_adj)
    if x_adj(i) >= 0
        C11 = 1;
        C21 = 0;
        if a_plus == 0
            a_plus=x_adj(i);
        end
    else
        C11 = 0;
        C21 = 1;
        if a_minus == 0
            a_minus=x_adj(i);
        end
    end
    if x_adj(i+1) >= 0
        C12 = 1;
        C22 = 0;
    else
        C12 = 0;
        C22 = 1;
    end
    if C11 == 1 && C12 ==1
        a_plus = a_plus +x_adj(i+1);
    elseif C21 == 1 && C22 ==1
        a_minus= a_minus + x_adj(i+1);
    elseif C11 == 1 && C12 ==0
        vec=[vec a_plus];
        a_plus=0;
        a_minus =0;
    elseif C21 == 1 && C22 ==0
        vec=[vec a_minus];
        a_plus=0;
        a_minus =0;
    end
    i=i+1;
end

dod_half= abs(vec);
eq_cyc_dod_100= zeros(N,1);
for i=1:length(vec)
    eq_cyc_dod_100(i,1) = 0.5*(dod_half(i))^(1.1);
end

total_cycle_100 = sum(eq_cyc_dod_100);