clear 
close all
clc

h=0.25;
network_name = 'threeBus_LV_industrial.m';

load('Building_Load.mat');
N = 96;
num_bus = 4;

load('day_ahead_calculated_battery_states.mat')
% load('case1_no_network_awareness.mat')
load('real_time_build_load.mat')
load('real_time_pv_gen.mat')
load('electrolyzer_load_active_reactive.mat')

pf_load = tan(acos(0.95));

P2 = 1*real_time_build_load;
Q2 = pf_load*P2;

P3 = 1*pv_gen;
Q3 = pf_load*P3*0;          % pv gen at unity power factor

%%

P_load = [P1 P2 -1*P3];
Q_load = [Q1 Q2 1*Q3];

% power factor
S_total = sqrt(sum(P_load').^2 + sum(Q_load').^2);
pf_load_agg_case0 = transpose(sum(P_load')./ S_total);

for time = 1:N
    define_constants;
    mpc = loadcase(network_name);
    mpc.bus(2,PD) = P_load(time,1);
    mpc.bus(3,PD) = P_load(time,2);
    mpc.bus(4,PD) = P_load(time,3);
    mpc.bus(2,QD) = Q_load(time,1);
    mpc.bus(3,QD) = Q_load(time,2);
    mpc.bus(4,QD) = Q_load(time,3);
    
    result_dg= runpf(mpc);
    Voltage_case0(:,time) = result_dg.bus(:,8);  
end

figure;plot(Voltage_case0')
hold on; yline(1.05); hold on; yline(0.95)

%% emission cost KPI

KPI5_case0 = 8.7622e+03;

KPI6_mean_abs_power_factor = mean(abs(pf_load_agg_case0));

KPI6_instances_below_lim = sum(abs(pf_load_agg_case0)<0.8)/(N)

%%

V_max = 1.05;
V_min = 0.95;
del_perm = 0.035;

s_agg1_post = sum(sum((Voltage_case0 > 1+del_perm).*(Voltage_case0 - 1 - del_perm)));
s_agg2_post = sum(sum((Voltage_case0 < 1-del_perm).*(1 - del_perm - Voltage_case0)));

l1=sum(sum(Voltage_case0 > V_max));
l2=sum(sum(Voltage_case0 > 1+del_perm));
l3=sum(sum(Voltage_case0 < 1-del_perm));
l4=sum(sum(Voltage_case0 < V_min));

%%

Voltage_node4 = Voltage_case0(4,:);

l1_n4=sum(sum(Voltage_node4 > V_max));
l2_n4=sum(sum(Voltage_node4 > 1+del_perm));
l3_n4=sum(sum(Voltage_node4 < 1-del_perm));
l4_n4=sum(sum(Voltage_node4 < V_min));

KPI4_voltage_PCC_Node4 = [l1_n4, l2_n4, l3_n4, l4_n4]/60

%% results
KPI1_profit_only_arbitrage = 0
KPI2_cycles_operation = 0

KPI3_CVC = [s_agg1_post, s_agg2_post]
KPI4_voltage_correction_index = [l1, l2, l3, l4]/4

KPI5_emission_cost = 8.7622e+03
KPI5_emission_saving = 0

KPI6_mean_abs_power_factor
KPI6_instances_below_lim

KPI4_voltage_PCC_Node4

%%

P_load_total = sum(P_load');
Q_load_total = sum(Q_load');
S_total_cal = sqrt((P_load_total).^2 + (Q_load_total).^2);

KPI7 = [max(P_load_total), min(P_load_total), var(P_load_total);...
    max(Q_load_total), min(Q_load_total), var(Q_load_total);...
    max(S_total_cal), min(S_total_cal), var(S_total_cal)]
    
KPI7_2 = [max(P_load_total), min(P_load_total), var(P_load_total),...
    max(Q_load_total), min(Q_load_total), var(Q_load_total),...
    max(S_total_cal), min(S_total_cal), var(S_total_cal)]';
%%

figure; 
subplot(211)
plot(0.25:0.25:24,P_load); hold on; plot(0.25:0.25:24,sum(P_load'))
subplot(212); 
plot(0.25:0.25:24,Q_load); hold on; plot(0.25:0.25:24,sum(Q_load'))

%%

C0_pf_load_agg = pf_load_agg_case0;
C0_Voltage = Voltage_case0;
% C0_Out_b_mat = Out_b_mat;
C0_P_load_total = sum(P_load')';
C0_Q_load_total = sum(Q_load')';