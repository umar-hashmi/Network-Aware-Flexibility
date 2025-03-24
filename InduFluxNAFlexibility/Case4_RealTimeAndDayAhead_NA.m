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

J_temp = -1*P3 + S_bat_case2(:,1);
P_load = [P1 P2 J_temp];
Q_load = [Q1 Q2 1*S_bat_case2(:,2)];

% power factor
S_total = sqrt(sum(P_load').^2 + sum(Q_load').^2);
pf_load_agg_post = transpose(sum(P_load')./ S_total);


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

P_T = sum(P_load')';
Q_T = sum(Q_load')';

figure; subplot(211); plot(P_T); hold on; plot(Q_T); 
subplot(212); plot(S_bat_case2(:,1)); hold on; plot(P3); hold on; plot(S_bat_case2(:,2))

figure; plot(pf_load_agg_post, '*')

%%

% battery parameters
e_ch=0.95;                          % Charging efficiency
e_dis =0.95;                        % Discharging efficiency 
del_max = 0.5/0.95;                     % Maximum charging rate
del_min = -0.5*0.95;                 % Minimum discharging rate
b_0 = 0.5;                         % Initial battery capacity
b_max = 1;                       % Maximum battery capacity
b_min = 0.2;                        % Minimum permissible battery capacity
h=0.25;                             % Sampling time

Q_b_max = 0.6;
S_b_max = 0.6;
pf_lim  = 0.8;

bat_para= [e_ch, e_dis, del_max del_min, b_max, b_min, S_b_max];

Delta_lim_P = 1*del_max;
Delta_lim_Q = 1*S_b_max;

J_temp = -1*P3 + S_bat_case2(:,1);
P_load = [P1 P2 J_temp];
Q_load = [Q1 Q2 1*S_bat_case2(:,2)];

V_max = 1.05;
V_min = 0.95;
del_perm = 0.035;

node_flexibility = 4;

index = 1;

for time = 1:N
    
    P_bat_da = S_bat_case2(time,1);
    Q_bat_da = S_bat_case2(time,2);
    P_L_rt = [P1(time,1) P2(time,1) -1*P3(time,1)+P_bat_da];
    Q_L_tr = [Q1(time,1) Q2(time,1) Q3(time,1)+Q_bat_da];
    
    del_b_old = 0;
    
    P_past = P_bat_da;
    Q_past = Q_bat_da;
    
    if time == 1
        b_past = b_0;
    end
    
    define_constants;
    mpc = loadcase(network_name);
    mpc.bus(2,PD) = P_L_rt(1,1);     mpc.bus(3,PD) = P_L_rt(1,2);     mpc.bus(4,PD) = P_L_rt(1,3);
    mpc.bus(2,QD) = Q_L_tr(1,1);     mpc.bus(3,QD) = Q_L_tr(1,2);     mpc.bus(4,QD) = Q_L_tr(1,3);
    result_dg= runpf(mpc);
    Voltage_temp = result_dg.bus(:,8);
    for t_inner = 1:15
        
        sampling_inner_loop = 1/60;
        
        % step 1: perform power flow
        if t_inner == 1
            [R_P_min4, R_P_max4, R_Q_min4, R_Q_max4] = Range_flex_operation(Voltage_temp(node_flexibility), 3, V_max, V_min, del_perm);
            [R_P_min3, R_P_max3, R_Q_min3, R_Q_max3] = Range_flex_operation(Voltage_temp(node_flexibility-1), 3, V_max, V_min, del_perm);
            R_P_min = max(R_P_min4, R_P_min3);
            R_P_max = min(R_P_max4,R_P_max3);
            R_Q_min = max(R_Q_min4, R_Q_min3);
            R_Q_max = min(R_Q_max4,R_Q_max3);
        else
            [R_P_min4, R_P_max4, R_Q_min4, R_Q_max4] = Range_flex_operation(Voltage_rt(node_flexibility,index-1), 3, V_max, V_min, del_perm);
            [R_P_min3, R_P_max3, R_Q_min3, R_Q_max3] = Range_flex_operation(Voltage_rt(node_flexibility-1,index-1), 3, V_max, V_min, del_perm);
            R_P_min = max(R_P_min4, R_P_min3);
            R_P_max = min(R_P_max4,R_P_max3);
            R_Q_min = max(R_Q_min4, R_Q_min3);
            R_Q_max = min(R_Q_max4,R_Q_max3);
        end
        
        % step 2: calculate real-time DOE
        R_P_min_corr = del_max*R_P_min;             R_P_max_corr = del_max*R_P_max;
        S_b_avail = sqrt(S_b_max^2 - (P_past).^2);
        Lim_inv = min(S_b_avail,Delta_lim_Q);
        if -1*S_b_avail*R_Q_min >= -1*S_b_avail*R_Q_max
            R_Q_max_corr = -1*S_b_avail*R_Q_min; R_Q_min_corr = -1*S_b_avail*R_Q_max;
        else
            R_Q_min_corr = -1*S_b_avail*R_Q_min;        R_Q_max_corr = -1*S_b_avail*R_Q_max;
        end
        
        % step 3: identify battery ranges
        ramp_slack_lower = (b_min - b_past)/sampling_inner_loop;     ramp_slack_upper = (b_max - b_past)/sampling_inner_loop;
        P_range_min = max(ramp_slack_lower, R_P_min_corr);           P_range_max = min(ramp_slack_upper, R_P_max_corr);
        
        % step 4: identify updated states of the battery
        [P_B_out, Q_B_out, del_bias, bat_ch_level] = LP_innerLoop(b_past, P_past, Q_past, del_b_old, P_range_min, P_range_max, R_Q_min_corr, R_Q_max_corr, bat_para);

        P_past = P_B_out;
        Q_past = Q_B_out;
        b_past = bat_ch_level;
        del_b_old = del_bias;
        
        Out_b_mat(index,:) = [P_B_out, Q_B_out, b_past, del_b_old];
        Envelope(index,:) = [P_range_min, P_range_max, R_Q_min_corr, R_Q_max_corr];
        
        define_constants;
        mpc = loadcase(network_name);
        mpc.bus(2,PD) = P_L_rt(1,1);     mpc.bus(3,PD) = P_L_rt(1,2);     mpc.bus(4,PD) = -1*P3(time,1)+P_B_out;
        mpc.bus(2,QD) = Q_L_tr(1,1);     mpc.bus(3,QD) = Q_L_tr(1,2);     mpc.bus(4,QD) = Q3(time,1)+Q_B_out;
        result_dg= runpf(mpc);
        Voltage_rt(:,index) = result_dg.bus(:,8);
        
        P_load_total(index,1) = P_L_rt(1,1)+P_L_rt(1,2)+-1*P3(time,1)+P_B_out;
        Q_load_total(index,1) = Q_L_tr(1,1)+Q_L_tr(1,2)+-1*Q3(time,1)+Q_B_out;
        % set b_past
        index = index +1;
    
    end
    
end
%%
t = 1/4:1/4:24;
t2 = 1/60:1/60:24;

figure; subplot(211); plot(t2,Envelope(:,1)); hold on; plot(t2,Envelope(:,2)); hold on; plot(t,S_bat_case2(:,1))
subplot(212); plot(t2,Envelope(:,3)); hold on; plot(t2,Envelope(:,4));
%%

figure
subplot(211)
plot(t,S_bat_case2(:,1))
hold on
plot(t2,Out_b_mat(:,1))
subplot(212)
plot(t,S_bat_case2(:,2))
hold on
plot(t2,Out_b_mat(:,2))
%%
load('case1_no_network_awareness.mat')

figure
subplot(211)
plot(t,S_bat_case1(:,1))
hold on
plot(t,S_bat_case2(:,1))
hold on
plot(t2,Out_b_mat(:,1))

subplot(212)
plot(t,S_bat_case1(:,2))
hold on
plot(t,S_bat_case2(:,2))
hold on
plot(t2,Out_b_mat(:,2))
%%

figure
curve1 = (1-del_perm)*ones(96,1);
curve2 = (1+del_perm)*ones(96,1);
hold on; shade(t,curve1,t,curve2,'FillType',[1 2;2 1])
hold on
curve1 = (1-del_perm)*ones(96,1);
curve2 = (V_min)*ones(96,1);
hold on; shade(t,curve1,t,curve2,'FillType',[1 2;2 1])
hold on
curve1 = (V_max)*ones(96,1);
curve2 = (1+del_perm)*ones(96,1);
hold on; shade(t,curve1,t,curve2,'FillType',[1 2;2 1])
hold on
curve1 = (V_max)*ones(96,1);
curve2 = (1.08)*ones(96,1);
hold on; shade(t,curve1,t,curve2,'FillType',[1 2;2 1])
hold on
curve1 = (0.92)*ones(96,1);
curve2 = (V_min)*ones(96,1);
hold on; shade(t,curve1,t,curve2,'FillType',[1 2;2 1])
hold on
plot(t,Voltage_post(4,:)')
hold on
plot(t2,Voltage_rt(4,:)')
hold on; yline(1.05); hold on; yline(0.95)
hold on; yline(1.035); hold on; yline(0.965)
title('voltage at node 4')
hold on




%%

figure
plot(t,Voltage_post(3,:)')
hold on
plot(t2,Voltage_rt(3,:)')
hold on; yline(1.05); hold on; yline(0.95)
hold on; yline(1.035); hold on; yline(0.965)
title('voltage at node 3')

%%

figure
subplot(211)
plot(t,Voltage_post')
hold on; yline(1.05); hold on; yline(0.95)
hold on; yline(1.035); hold on; yline(0.965)
subplot(212)
plot(t2,Voltage_rt')
hold on; yline(1.05); hold on; yline(0.95)
hold on; yline(1.035); hold on; yline(0.965)

%%

figure
stairs([0,t],[Voltage_post(4,1),Voltage_post(4,:)])
hold on
stairs([0,t2],[Voltage_rt(4,1),Voltage_rt(4,:)])
hold on; yline(1.05); hold on; yline(0.95)
hold on; yline(1.035); hold on; yline(0.965)
%% Profit calculation

load('data_nn.mat')
kappa = 1;

carbon_intensity=Dataneeded(:,1); %emission cost g/kWh
day_ahead_price=Dataneeded(:,2); %DAM cost euros/MWh
emission_cost = 100; %euros/tonne
DA_price_1_day = resample(day_ahead_price(1:24), 4,1);
CI_cost_1_day = resample(carbon_intensity(1:24), 4,1);
DA_price_1_day_1min = resample(day_ahead_price(1:24), 60,1);
CI_cost_1_day_1min = resample(carbon_intensity(1:24), 60,1);
Total_cost_DA = DA_price_1_day + CI_cost_1_day*emission_cost/10e3;
Total_cost_DA_1min = DA_price_1_day_1min + CI_cost_1_day_1min*emission_cost/10e3;
% figure; plot(Total_cost_DA); hold on; plot(DA_price_1_day); hold on; plot(CI_cost_1_day*emission_cost/10e3)


profit_only_arbitrage_dayAhead =  0.25*(sum(Total_cost_DA'*subplus(S_bat_case2(:,1))-kappa*Total_cost_DA'*subplus(-S_bat_case2(:,1))))
profit_only_arbitrage_RT =  (sum(Total_cost_DA_1min'*subplus(Out_b_mat(:,1))-kappa*Total_cost_DA_1min'*subplus(-Out_b_mat(:,1))*e_dis))/60

%% emission cost KPI
% P_total = sum(P_load')';


KPI5_emission_cost = sum(CI_cost_1_day_1min.*P_load_total)/60;

KPI5_case0 = 8.7622e+03;

KPI5_emission_saving = KPI5_case0 - KPI5_emission_cost;

%%

pf_load_agg_case3 = transpose((P_load_total)./ sqrt((P_load_total).^2 + (Q_load_total).^2));

KPI6_mean_abs_power_factor = mean(abs(pf_load_agg_case3));

KPI6_instances_below_lim = sum(abs(pf_load_agg_case3)<0.8)/(24*60)

figure; subplot(211); plot(t2,P_load_total); hold on; plot(t2,Q_load_total); subplot(212); plot(t2,pf_load_agg_case3)
%%
x_out_bat_case3 = diff(Out_b_mat(:,3));

s_agg1_post = sum(sum((Voltage_rt > 1+del_perm).*(Voltage_rt - 1 - del_perm)));
s_agg2_post = sum(sum((Voltage_rt < 1-del_perm).*(1 - del_perm - Voltage_rt)));


l1=sum(sum(Voltage_rt > V_max));
l2=sum(sum(Voltage_rt > 1+del_perm));
l3=sum(sum(Voltage_rt < 1-del_perm));
l4=sum(sum(Voltage_rt < V_min));

%%

Voltage_node4 = Voltage_rt(4,:);

l1_n4=sum(sum(Voltage_node4 > V_max));
l2_n4=sum(sum(Voltage_node4 > 1+del_perm));
l3_n4=sum(sum(Voltage_node4 < 1-del_perm));
l4_n4=sum(sum(Voltage_node4 < V_min));

KPI4_voltage_PCC_Node4 = [l1_n4, l2_n4, l3_n4, l4_n4]/60

%%

% results
KPI1_profit_only_arbitrage = profit_only_arbitrage_RT
KPI2_cycles_operation = calculate_cycles(x_out_bat_case3)

KPI3_CVC = [s_agg1_post, s_agg2_post]
KPI4_voltage_correction_index = [l1, l2, l3, l4]/60

KPI5_emission_cost
KPI5_emission_saving

KPI6_mean_abs_power_factor
KPI6_instances_below_lim

%%

S_total_cal = sqrt((P_load_total).^2 + (Q_load_total).^2);

KPI7 = [max(P_load_total), min(P_load_total), var(P_load_total);...
    max(Q_load_total), min(Q_load_total), var(Q_load_total);...
    max(S_total_cal), min(S_total_cal), var(S_total_cal)]
    
KPI7_2 = [max(P_load_total), min(P_load_total), var(P_load_total),...
    max(Q_load_total), min(Q_load_total), var(Q_load_total),...
    max(S_total_cal), min(S_total_cal), var(S_total_cal)]';

%%

C4_pf_load_agg = pf_load_agg_case3;
C4_Voltage = Voltage_rt;
C4_Out_b_mat = Out_b_mat;
C4_P_load_total = P_load_total;
C4_Q_load_total = Q_load_total;
C4_battery_ch_level = Out_b_mat(:,3);
% C4_Q_load_total = Q_load_total;

%%

Voltage_rt_resample_N4 = downsample(Voltage_rt(4,:)',15);
Voltage_rt_resample_N3 = downsample(Voltage_rt(3,:)',15);

Pflex = downsample(Out_b_mat(:,1),15);
Qflex = downsample(Out_b_mat(:,2),15);
