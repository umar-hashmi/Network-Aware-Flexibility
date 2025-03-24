clear 
close all
clc

%% Utilizing flexibility for committments in day ahead market
t = 0.25:0.25:24;
% battery parameters
e_ch=0.95;                          % Charging efficiency
e_dis =0.95;                        % Discharging efficiency 
del_max = 0.5;                     % Maximum charging rate
del_min = -del_max;                 % Minimum discharging rate
b_0 = 0.5;                         % Initial battery capacity
b_max = 1;                       % Maximum battery capacity
b_min = 0.2;                        % Minimum permissible battery capacity
h=0.25;                             % Sampling time

Q_b_max = 0.6;
S_b_max = 0.6;
pf_lim  = 0.8;

load('data_nn.mat')

carbon_intensity=Dataneeded(:,1); %emission cost g/kWh
day_ahead_price=Dataneeded(:,2); %DAM cost euros/MWh
emission_cost = 100; %euros/tonne

%%

DA_price_1_day = resample(day_ahead_price(1:24), 4,1);
CI_cost_1_day = resample(carbon_intensity(1:24), 4,1);
Total_cost_DA = DA_price_1_day + CI_cost_1_day*emission_cost/10e3;

figure; plot(Total_cost_DA); hold on; plot(DA_price_1_day); hold on; plot(CI_cost_1_day*emission_cost/10e3)

%%

kappa = 1;                          % the ratio of selling price and buying price   
N=96;

real_buy = Total_cost_DA/e_ch;
real_sell = kappa*Total_cost_DA*e_dis;

A_buy = diag(real_buy);
A_sell = diag(real_sell);
A_minus = -1*eye(N);
A_zero = zeros(N);

%%
tic
% LP matrix formulation
A = [A_buy A_minus; A_sell A_minus; tril(ones(N,N),-1) + eye(N)  A_zero; -tril(ones(N,N),-1) - eye(N)  A_zero];
b=[zeros(2*N,1);(b_max-b_0)*ones(N,1); (b_0-b_min)*ones(N,1)];

lb=[del_min*ones(N,1)*h; -100000*ones(N,1)];
ub=[del_max*ones(N,1)*h; 100000*ones(N,1)];

Aeq=[];
beq=[];
f=[zeros(N,1); ones(N,1)];

x_state = linprog(f,A,b,Aeq,beq,lb,ub);
toc
battery_ramp= x_state(1:N);

KPI1_profit_only_arbitrage =  sum(Total_cost_DA'*subplus(battery_ramp)/e_ch-kappa*Total_cost_DA'*subplus(-battery_ramp)*e_dis);

figure; plot(t,battery_ramp)

x_adj= battery_ramp/b_max;

%%

x_ch = max(0,battery_ramp);
x_ds = -min(0,battery_ramp);
x_bat_out =x_ch/e_ch - x_ds*e_dis;
bat_cap = b_0+cumsum(battery_ramp);
figure; plot(t, bat_cap)

P_bat = 4*x_bat_out;
Q_bat = zeros(N,1);
%%

S_bat_case1 = [P_bat, Q_bat];
figure; plot(t, S_bat_case1)
save 'case1_no_network_awareness.mat' S_bat_case1

%%

network_name = 'threeBus_LV_industrial.m';

load('Building_Load.mat');
N = 96;
num_bus = 4;

% load('day_ahead_calculated_battery_states.mat')
load('real_time_build_load.mat')
load('real_time_pv_gen.mat')
load('electrolyzer_load_active_reactive.mat')

pf_load = tan(acos(0.95));
P2 = 1*real_time_build_load;
Q2 = pf_load*P2;

P3 = 1*pv_gen;
Q3 = pf_load*P3*0;          % pv gen at unity power factor

J_temp = -1*P3 + S_bat_case1(:,1);
P_load = [P1 P2 J_temp];
Q_load = [Q1 Q2 1*S_bat_case1(:,2)];

% power factor
S_total = sqrt(sum(P_load').^2 + sum(Q_load').^2);
pf_load_agg_case1 = transpose(sum(P_load')./ S_total);

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
    Voltage_post(:,time) = result_dg.bus(:,8); 
end

figure;plot(Voltage_post')
hold on; yline(1.05); hold on; yline(0.95)

%%

figure; plot(pf_load_agg_case1, '*')

KPI6_mean_abs_power_factor = mean(abs(pf_load_agg_case1));

KPI6_instances_below_lim = sum(abs(pf_load_agg_case1)<0.8)/N

%%

KPI2_cycles_operation = calculate_cycles(x_bat_out);

%% voltage correction metric

V_max = 1.05;
V_min = 0.95;
del_perm = 0.035;

s_agg1_post = sum(sum((Voltage_post > 1+del_perm).*(Voltage_post - 1 - del_perm)));
s_agg2_post = sum(sum((Voltage_post < 1-del_perm).*(1 - del_perm - Voltage_post)));

KPI3_CVC = [s_agg1_post, s_agg2_post]

l1=sum(sum(Voltage_post > V_max));
l2=sum(sum(Voltage_post > 1+del_perm));
l3=sum(sum(Voltage_post < 1-del_perm));
l4=sum(sum(Voltage_post < V_min));

KPI4_voltage_correction_index = 0.25*[l1, l2, l3, l4]

%%

Voltage_node4 = Voltage_post(4,:);

l1_n4=sum(sum(Voltage_node4 > V_max));
l2_n4=sum(sum(Voltage_node4 > 1+del_perm));
l3_n4=sum(sum(Voltage_node4 < 1-del_perm));
l4_n4=sum(sum(Voltage_node4 < V_min));

KPI4_voltage_PCC_Node4 = [l1_n4, l2_n4, l3_n4, l4_n4]/60

%% emission cost KPI
P_total = sum(P_load')';

KPI5_emission_cost = 0.25*sum(CI_cost_1_day.*P_total);
KPI5_case0 = 0.25*sum(CI_cost_1_day.*sum([P1 P2 -P3]')');
KPI5_emission_saving = KPI5_case0 - KPI5_emission_cost;

%% results
KPI1_profit_only_arbitrage
KPI2_cycles_operation
KPI3_CVC
KPI4_voltage_correction_index

KPI5_emission_cost
KPI5_emission_saving

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

C1_pf_load_agg = pf_load_agg_case1;
C1_Voltage = Voltage_post;
C1_Out_b_mat = S_bat_case1;
C1_P_load_total = sum(P_load')';
C1_Q_load_total = sum(Q_load')';
C1_battery_ch_level = bat_cap;
