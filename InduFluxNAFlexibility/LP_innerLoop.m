function [P_B_out, Q_B_out, del_bias, bat_ch_level] = LP_innerLoop(b_past, P_bat, Q_bat, del_b_old, Pmin, Pmax, Qmin, Qmax, bat_para)

% input battery capacity, renewable generation
e_ch=bat_para(1);                          % Charging efficiency
e_dis =bat_para(2);                        % Discharging efficiency 
del_max = bat_para(3)/e_ch;                     % Maximum charging rate
del_min = bat_para(4)*e_dis;                 % Minimum discharging rate
b_max = bat_para(5);                       % Maximum battery capacity
b_min = bat_para(6);                        % Minimum permissible battery capacity
S_max = bat_para(7);

k= 1/60;

%% Active power 

A=[1, -1; -1, -1];
% b=[P_bat - del_b_old; -P_bat + del_b_old];
b=[P_bat; -P_bat];
f=[0; 1];

Aeq=[];
beq=[];
lb = [Pmin ; -1000000];
ub = [Pmax ; 1000000];

[x_state,fval,exitflag,output] = linprog(f,A,b,Aeq,beq,lb,ub);


if exitflag == 1
    P_B_out = x_state(1);
    del_bias =  x_state(2);
else
    P_B_out = 0;
    del_bias =  del_b_old;
end

x_ch = max(0,P_B_out);
x_ds = -min(0,P_B_out);
x_bat_out =x_ch*e_ch - x_ds/e_dis;
bat_ch_level = b_past+x_bat_out*k;


%% Reactive power

A=[1, -1; -1, -1];
b=[Q_bat ; -Q_bat ];
f=[0; 1];

Aeq=[];
beq=[];
lb= [Qmin; -1000000];
ub= [Qmax; 1000000];

[x_state,fval,exitflag,output] = linprog(f,A,b,Aeq,beq,lb,ub);

if exitflag == 1
    Q_B_out = x_state(1);
else
    Q_B_out = 0;
end


