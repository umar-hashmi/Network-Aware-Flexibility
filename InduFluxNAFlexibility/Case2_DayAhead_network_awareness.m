clear
close all
clc

h=0.25;
network_name = 'threeBus_LV_industrial.m';
load('Building_Load.mat');

N = 96;
num_bus = 4;

%% 362 days of load 15 min data

build_load = resultsbuildwbatt(:,2);
day_ahead_price = resultsbuildwbatt(:,3);

buildLoad_reshape= 1.5*reshape(build_load, 96,362)/1000;        %changed from 1.8

%% electrolyzer load profile

load('Active_Power_electrolyzer_V2.mat')
load('Reactive_Power_electrolyzer_V2.mat')
P2=load('Var_constant_V1.mat');

t = 0.25:0.25:24;
tsin = timeseries(Q.Data,Q.Time);
tsin_Q = timeseries(Q1.Data,Q1.Time);
tsout_P = resample(tsin,t);
tsout_Q = resample(tsin_Q,t);
tsOUT_P = tsout_P.Data;
tsOUT_Q = tsout_Q.Data;
tsin_P2 = timeseries(P2.Q.Data,P2.Q.Time);
tsout_P2 = resample(tsin_P2,t);
tsOUT_P2 = tsout_P2.Data;
figure; plot(t, 2*tsout_P.Data./max(tsOUT_P)); hold on; plot(t,8*tsout_Q.Data./max(tsOUT_P))    %;hold on; plot(tsOUT_P2);

%% PV generation profile

pv_size = 1.5; %in MWp  changed from 2 to 1.5

normalized_solar = csvread('SolarPV_normalized_profiles.csv');
tsin1 = timeseries(normalized_solar,1:24);
normalized_solar_resampled = resample(tsin1,t');
n_s_data = normalized_solar_resampled.Data;
n_s_data(isnan(n_s_data))=0;
day_ahead_solar = pv_size*n_s_data;

%%
pf_load = tan(acos(0.95));

P1 = 2*tsOUT_P./max(tsOUT_P);           % changed from 2.5 to 2
Q1 = 8*tsOUT_Q./max(tsOUT_P);           % changed from 10 to 8

P2 = buildLoad_reshape;
Q2 = pf_load*P2;

P3 = day_ahead_solar;
Q3 = pf_load*P3*0;  % pv gen at unity power factor

Voltage = zeros(100,4,96);

sen_max= 100;

for scenario =1:sen_max
    
    P_load = [P1 P2(:,scenario) -P3(:,scenario)];
    Q_load = [Q1 Q2(:,scenario) Q3(:,scenario)];
    
    % power factor
    S_total = sqrt(sum(P_load').^2 + sum(Q_load').^2);
    pf_load_agg(:,scenario) = transpose(sum(P_load')./ S_total);
    P_total_uncorrected(:,scenario) = sum(P_load')';
    Q_total_uncorrected(:,scenario) = sum(Q_load')';
    
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
        Voltage(scenario,:,time) = result_dg.bus(:,8);
        
    end
end

%%
V_low = zeros(4,96);
V_up = zeros(4,96);
for i = 1:4
    for j = 1:96
        
        V_new = squeeze(Voltage(:,i,j));
        V_low(i,j) = ecdf_quartile(V_new,0.05);
        V_up(i,j) = ecdf_quartile(V_new,0.95);

    end
end

%%
V_max = 1.05;
V_min = 0.95;
del_perm = 0.035;

P_range_ANRC_lower = zeros(4,96);
P_range_ANRC_upper = zeros(4,96);
Q_range_PRC_lower = zeros(4,96);
Q_range_PRC_upper = zeros(4,96);

for i = 1:4
    for  j =1:96
        [R_P_min_l, R_P_max_l, R_Q_min_l, R_Q_max_l] = Range_flex_operation(V_low(i,j), 3, V_max, V_min, del_perm);
        [R_P_min_u, R_P_max_u, R_Q_min_u, R_Q_max_u] = Range_flex_operation(V_up(i,j), 3, V_max, V_min, del_perm);

        P_range_ANRC_lower(i,j)= max(R_P_min_l,R_P_min_u);
        P_range_ANRC_upper(i,j)= min(R_P_max_l,R_P_max_u);

        Q_range_PRC_lower(i,j)= max(R_Q_min_l,R_Q_min_u);
        Q_range_PRC_upper(i,j)= min(R_Q_max_l,R_Q_max_u);
    end
end

%%

t=t';
node = 4;

figure
subplot(311)
curve1 = (1-del_perm)*ones(96,1);
curve2 = (1+del_perm)*ones(96,1);
hold on; shade(t,curve1,t,curve2,'FillType',[1 2;2 1])
hold on
plot(t, V_low(node,:)); hold on; plot(t, V_max*ones(96,1));hold on; plot(t, V_min*ones(96,1));
hold on; plot(t, V_up(node,:))
subplot(312)
curve1 = P_range_ANRC_lower(node,:);
curve2 = P_range_ANRC_upper(node,:);
shade(t,curve1,t,curve2,'FillType',[1 2;2 1],'LineWidth',2);
subplot(313)
curve1 = Q_range_PRC_lower(node,:);
curve2 = Q_range_PRC_upper(node,:);
shade(t,curve1,t,curve2,'FillType',[1 2;2 1],'LineWidth',2);


%% Utilizing flexibility for committments in day ahead market

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

real_buy = Total_cost_DA/e_ch;
real_sell = kappa*Total_cost_DA*e_dis;

A_buy = diag(real_buy);
A_sell = diag(real_sell);
A_minus = -1*eye(N);
A_zero = zeros(N);

%%
node = 4;
tic
% LP matrix formulation
A = [A_buy A_minus; A_sell A_minus; tril(ones(N,N),-1) + eye(N)  A_zero; -tril(ones(N,N),-1) - eye(N)  A_zero];
b=[zeros(2*N,1);(b_max-b_0)*ones(N,1); (b_0-b_min)*ones(N,1)];

lb=[max(del_min*ones(N,1), del_max*transpose(P_range_ANRC_lower(node,:)))*h; -100000*ones(N,1)];
ub=[min(del_max*ones(N,1), del_max*transpose(P_range_ANRC_upper(node,:)))*h; 100000*ones(N,1)];

Aeq=[];
beq=[];
f=[zeros(N,1); ones(N,1)];

x_state = linprog(f,A,b,Aeq,beq,lb,ub);
toc
battery_ramp= x_state(1:N);

KPI1_profit_only_arbitrage =  sum(Total_cost_DA'*subplus(battery_ramp)/e_ch-kappa*Total_cost_DA'*subplus(-battery_ramp)*e_dis)

figure; plot(t,battery_ramp)

x_adj= battery_ramp/b_max;

%%

x_ch = max(0,battery_ramp);
x_ds = -min(0,battery_ramp);
x_bat_out =x_ch/e_ch - x_ds*e_dis;
bat_cap = b_0+cumsum(battery_ramp);
figure; plot(t, bat_cap)

%% reactive power setpoints TO IMPROVE THIS
% consider the power factor too here

%
figure
boxplot(pf_load_agg')
hold on 
plot(mean(pf_load_agg'),'*')

pf_mean_uncorrected = mean(pf_load_agg')';
Q_T_mean = mean(Q_total_uncorrected')';
P_T_mean = mean(P_total_uncorrected')';

figure; plot(Q_T_mean); hold on; plot(P_T_mean)
%%
S_b_avail = sqrt(S_b_max^2 - (abs(x_bat_out)*4).^2);

Q_l_lim = Q_range_PRC_lower*Q_b_max;
Q_u_lim = Q_range_PRC_upper*Q_b_max;

for i = 1:96
    if Q_range_PRC_lower(node,i) > 0
        Q_bat(i,1) = -1*min(S_b_avail(i), Q_range_PRC_lower(node,i));
    elseif Q_range_PRC_upper(node,i) < 0
        Q_bat(i,1) = -1*max(-S_b_avail(i), Q_range_PRC_upper(node,i));
    else       
        if abs(pf_mean_uncorrected(i,1)) < 0.85
            if abs(P_T_mean(i,1)) > 0.75
                if Q_T_mean(i,1) >= 0
                    Q_bat(i,1) = min(max(-S_b_avail(i), 0.75*(abs(P_T_mean(i,1))) - Q_T_mean(i,1)), 0);
                else
                    Q_bat(i,1) = max(min(S_b_avail(i), -0.75*(abs(P_T_mean(i,1))) - Q_T_mean(i,1)), 0);
                end
            else
                if Q_load(i,1) >= 0
                    Q_bat(i,1) = min(max(-S_b_avail(i), - Q_T_mean(i,1)), 0);
                else
                    Q_bat(i,1) = max(min(S_b_avail(i),  - Q_T_mean(i,1)), 0);
                end
            end
        else 
            Q_bat(i,1) = 0;
        end
        
    end
end

figure; 
plot(Q_bat);hold on; plot(S_b_avail)

%% Veryfying the battery states

pf_load = tan(acos(0.95));
Voltage_post = zeros(100,4,96);

for scenario =1:sen_max
    J_temp = -1*P3(:,scenario)+ 4*x_bat_out;
    P_load = [P1 P2(:,scenario) J_temp];
    Q_load = [Q1 Q2(:,scenario) 1*Q_bat];   
    % power factor
    S_total = sqrt(sum(P_load').^2 + sum(Q_load').^2);
    pf_load_agg_post(:,scenario) = transpose(sum(P_load')./ S_total);
     
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
        Voltage_post(scenario,:,time) = result_dg.bus(:,8);        
    end
end

V_low = zeros(4,96);
V_up = zeros(4,96);
for i = 1:4
    for j = 1:96        
        V_new = squeeze(Voltage_post(:,i,j));
        V_low(i,j) = ecdf_quartile(V_new,0.05);
        V_up(i,j) = ecdf_quartile(V_new,0.95);
    end
end

P_range_ANRC_lower_post = zeros(4,96);
P_range_ANRC_upper_post = zeros(4,96);
Q_range_PRC_lower_post = zeros(4,96);
Q_range_PRC_upper_post = zeros(4,96);

for i = 1:4
    for  j =1:96
        [R_P_min_l, R_P_max_l, R_Q_min_l, R_Q_max_l] = Range_flex_operation(V_low(i,j), 3, V_max, V_min, del_perm);
        [R_P_min_u, R_P_max_u, R_Q_min_u, R_Q_max_u] = Range_flex_operation(V_up(i,j), 3, V_max, V_min, del_perm);

        P_range_ANRC_lower_post(i,j)= max(R_P_min_l,R_P_min_u);
        P_range_ANRC_upper_post(i,j)= min(R_P_max_l,R_P_max_u);

        Q_range_PRC_lower_post(i,j)= max(R_Q_min_l,R_Q_min_u);
        Q_range_PRC_upper_post(i,j)= min(R_Q_max_l,R_Q_max_u);
    end
end

t=t';

figure
subplot(311)
curve1 = (1-del_perm)*ones(96,1);
curve2 = (1+del_perm)*ones(96,1);
hold on; shade(t,curve1,t,curve2,'FillType',[1 2;2 1])
hold on
plot(t, V_low(node,:)); hold on; plot(t, V_max*ones(96,1));hold on; plot(t, V_min*ones(96,1));
hold on; plot(t, V_up(node,:))
subplot(312)
curve1 = P_range_ANRC_lower_post(node,:);
curve2 = P_range_ANRC_upper_post(node,:);
shade(t,curve1,t,curve2,'FillType',[1 2;2 1],'LineWidth',2);
subplot(313)
curve1 = Q_range_PRC_lower_post(node,:);
curve2 = Q_range_PRC_upper_post(node,:);
shade(t,curve1,t,curve2,'FillType',[1 2;2 1],'LineWidth',2);

%%
figure
plot(mean(abs(pf_load_agg'))')
hold on
plot(mean(abs(pf_load_agg_post'))','*-')
hold on
yline(0.8)
hold on
yline(-0.8)
%%
pf_mean_corrected = mean(pf_load_agg_post')';

figure
plot(pf_mean_uncorrected)
hold on
plot(pf_mean_corrected,'*-')
hold on
yline(0.8)
hold on
yline(-0.8)


%% gauging voltage envelope

R_P_post  = sum(P_range_ANRC_upper_post(node,:) - P_range_ANRC_lower_post(node,:));
R_P  = sum(P_range_ANRC_upper(node,:) - P_range_ANRC_lower(node,:));

Improvement_P = 100*(R_P_post-R_P)/R_P

R_Q_post  = sum(Q_range_PRC_upper_post(node,:) - Q_range_PRC_lower_post(node,:));
R_Q  = sum(Q_range_PRC_upper(node,:) - Q_range_PRC_lower(node,:));

Improvement_Q = 100*(R_Q_post-R_Q)/R_Q

%% Voltage correction index
node = 4;
s1=sum(sum(Voltage(:,node,:) > V_max));
s2=sum(sum(Voltage(:,node,:) > 1+del_perm));
s3=sum(sum(Voltage(:,node,:) < 1-del_perm));
s4=sum(sum(Voltage(:,node,:) < V_min));

s = [s1, s2, s3, s4]

l1=sum(sum(Voltage_post(:,node,:) > V_max));
l2=sum(sum(Voltage_post(:,node,:) > 1+del_perm));
l3=sum(sum(Voltage_post(:,node,:) < 1-del_perm));
l4=sum(sum(Voltage_post(:,node,:) < V_min));

l = [l1, l2, l3, l4]

%% Cumulative voltage correction

s_agg1 = sum(sum((Voltage(:,node,:) > 1+del_perm).*(Voltage(:,node,:) - 1 - del_perm)));
s_agg2 = sum(sum((Voltage(:,node,:) < 1-del_perm).*(1 - del_perm - Voltage(:,node,:))));

s_agg1_post = sum(sum((Voltage_post(:,node,:) > 1+del_perm).*(Voltage_post(:,node,:) - 1 - del_perm)));
s_agg2_post = sum(sum((Voltage_post(:,node,:) < 1-del_perm).*(1 - del_perm - Voltage_post(:,node,:))));

CVC = [s_agg1, s_agg2; s_agg1_post, s_agg2_post]

%%

figure
subplot(211)
histogram(pf_load_agg)
subplot(212)
histogram(pf_load_agg_post)

pf_compare = [mean(mean(pf_load_agg)) mean(mean(pf_load_agg_post))]

%%
P_T = zeros(96,sen_max);
Q_T = zeros(96,sen_max);
S_T = zeros(96,sen_max);
for scenario =1:sen_max
    J_temp = -1*P3(:,scenario)+ 4*x_bat_out;
    P_load = [P1 P2(:,scenario) J_temp];
    Q_load = [Q1 Q2(:,scenario) -1*Q_bat];
    
    % power factor
    S_total = sqrt(sum(P_load').^2 + sum(Q_load').^2);
    pf_load_agg_post(:,scenario) = transpose(sum(P_load')./ S_total);
    
    P_T(:,scenario) = sum(P_load')';
    Q_T(:,scenario) = sum(Q_load')';
    S_T(:,scenario) = S_total';
end

figure
subplot(311); plot(P_T)
subplot(312); plot(Q_T)
subplot(313); plot(S_T)
%%

P_bat_save = 4*x_bat_out;
Q_bat_save = -1*Q_bat;
S_bat_case2 = [P_bat_save, Q_bat_save];
figure; plot(t, S_bat_case2)

%%
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

J_temp = -1*P3 + S_bat_case2(:,1);
P_load = [P1 P2 J_temp];
Q_load = [Q1 Q2 1*S_bat_case2(:,2)];

% power factor
S_total = sqrt(sum(P_load').^2 + sum(Q_load').^2);
pf_load_agg_case2 = transpose(sum(P_load')./ S_total);

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
    Voltage_case2(:,time) = result_dg.bus(:,8); 
end

figure;plot(Voltage_case2')
hold on; yline(1.05); hold on; yline(0.95)

%%

V_max = 1.05;
V_min = 0.95;
del_perm = 0.035;

s_agg1_post = sum(sum((Voltage_case2 > 1+del_perm).*(Voltage_case2 - 1 - del_perm)));
s_agg2_post = sum(sum((Voltage_case2 < 1-del_perm).*(1 - del_perm - Voltage_case2)));

l1=sum(sum(Voltage_case2 > V_max));
l2=sum(sum(Voltage_case2 > 1+del_perm));
l3=sum(sum(Voltage_case2 < 1-del_perm));
l4=sum(sum(Voltage_case2 < V_min));

%%

Voltage_node4 = Voltage_case2(4,:);

l1_n4=sum(sum(Voltage_node4 > V_max));
l2_n4=sum(sum(Voltage_node4 > 1+del_perm));
l3_n4=sum(sum(Voltage_node4 < 1-del_perm));
l4_n4=sum(sum(Voltage_node4 < V_min));

KPI4_voltage_PCC_Node4 = [l1_n4, l2_n4, l3_n4, l4_n4]/60


%% emission cost KPI
P_total = sum(P_load')';

KPI5_emission_cost = sum(CI_cost_1_day.*P_total)*0.25;

KPI5_case0 = 0.25*sum(CI_cost_1_day.*sum([P1 P2 -P3]')');

KPI5_emission_saving = KPI5_case0 - KPI5_emission_cost;

%% Power factor KPI

KPI6_mean_abs_power_factor = mean(abs(pf_load_agg_case2));

KPI6_instances_below_lim = sum(abs(pf_load_agg_case2)<0.8)/N

%%

% results
KPI1_profit_only_arbitrage =  sum(Total_cost_DA'*subplus(battery_ramp)/e_ch-kappa*Total_cost_DA'*subplus(-battery_ramp)*e_dis)

KPI2_cycles_operation = calculate_cycles(x_bat_out)
KPI3_CVC = [s_agg1_post, s_agg2_post]
KPI4_voltage_correction_index = 0.25*[l1, l2, l3, l4]

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
x_ch = max(0,battery_ramp);
x_ds = -min(0,battery_ramp);
x_bat_out =x_ch/e_ch - x_ds*e_dis;
bat_cap = b_0+cumsum(battery_ramp);
figure; plot(t, bat_cap)
%%

C2_pf_load_agg = pf_load_agg_case2;
C2_Voltage = Voltage_case2;
C2_Out_b_mat = S_bat_case2;
C2_P_load_total = sum(P_load')';
C2_Q_load_total = sum(Q_load')';
C2_battery_ch_level = bat_cap;
