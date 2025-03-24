function [R_P_min, R_P_max, R_Q_min, R_Q_max] = Range_flex_operation(V_k, mode, V_max, V_min, del_perm)

P_max = 1;
P_min = -1;
Q_max = 1;
Q_min = -1;

%mode 1 =  PRC for P, Q
%mode 2 =  ANRC for P, Q
%mode 3 = ANRC for P and RPC for Q
%mode 4 = only Q control with PRC
%mode 5 = only Q control with ANRC
%mode otherwise no limit of regions

if mode == 1 % Positive reinforcement voltage regulation inverter op
    if V_k < V_min
        R_P_min = P_min;
        R_P_max = P_min;
        R_Q_min = Q_max;
        R_Q_max = Q_max;
    elseif V_k >= V_min && V_k < 1- del_perm
        R_P_min = P_min;
        R_P_max = P_min*(V_k - (1-del_perm))/(V_min - (1-del_perm));
        R_Q_min = Q_max*(V_k - (1-del_perm))/(V_min - (1-del_perm));
        R_Q_max = Q_max;     
    elseif V_k >= 1- del_perm && V_k <= 1+ del_perm
        R_P_min = P_min;
        R_P_max = P_max;
        R_Q_min = Q_min;
        R_Q_max = Q_max;
    elseif V_k > 1+ del_perm && V_k <= V_max
        R_P_min = P_max*(V_k - (1+del_perm))/(V_max - (1+del_perm));
        R_P_max = P_max;
        R_Q_min = Q_min;
        R_Q_max = Q_min*(V_k - (1+del_perm))/(V_max - (1+del_perm));
    elseif V_k > V_max
        R_P_min = P_max;
        R_P_max = P_max;
        R_Q_min = Q_min;
        R_Q_max = Q_min;
    end
    
elseif mode == 2
    if V_k < V_min
        R_P_min = P_min;
        R_P_max = 0;
        R_Q_min = 0;
        R_Q_max = Q_max;
    elseif V_k >= V_min && V_k < 1- del_perm
        R_P_min = P_min;
        R_P_max = P_max*(V_min - V_k )/(V_min - (1-del_perm));
        R_Q_min = Q_min*(V_min - V_k)/(V_min - (1-del_perm));
        R_Q_max = Q_max;     
    elseif V_k >= 1- del_perm && V_k <= 1+ del_perm
        R_P_min = P_min;
        R_P_max = P_max;
        R_Q_min = Q_min;
        R_Q_max = Q_max;
    elseif V_k > 1+ del_perm && V_k <= V_max
        R_P_min = P_min*(V_max - V_k )/(V_max - (1+del_perm));
        R_P_max = P_max;
        R_Q_min = Q_min;
        R_Q_max = Q_max*(V_max - V_k)/(V_max - (1+del_perm));
    elseif V_k > V_max
        R_P_min = 0;
        R_P_max = P_max;
        R_Q_min = Q_min;
        R_Q_max = 0;
    end
    
elseif mode == 3    % Hybrid P by ANRC, Q by PRC
    if V_k < V_min
        R_P_min = P_min;
        R_P_max = 0;
        R_Q_min = Q_max;
        R_Q_max = Q_max;
    elseif V_k >= V_min && V_k < 1- del_perm
        R_P_min = P_min;
        R_P_max = P_max*(V_min - V_k )/(V_min - (1-del_perm));
        R_Q_min = Q_max*(V_k - (1-del_perm))/(V_min - (1-del_perm));
        R_Q_max = Q_max;     
    elseif V_k >= 1- del_perm && V_k <= 1+ del_perm
        R_P_min = P_min;
        R_P_max = P_max;
        R_Q_min = Q_min;
        R_Q_max = Q_max;
    elseif V_k > 1+ del_perm && V_k <= V_max
        R_P_min = P_min*(V_max - V_k )/(V_max - (1+del_perm));
        R_P_max = P_max;
        R_Q_min = Q_min;
        R_Q_max = Q_min*(V_k - (1+del_perm))/(V_max - (1+del_perm));
    elseif V_k > V_max
        R_P_min = 0;
        R_P_max = P_max;
        R_Q_min = Q_min;
        R_Q_max = Q_min;
    end
    
elseif mode == 4        %Q(V) PRC
    if V_k < V_min
        R_P_min = P_min;
        R_P_max = P_max;
        R_Q_min = Q_max;
        R_Q_max = Q_max;
    elseif V_k >= V_min && V_k < 1- del_perm
        R_P_min = P_min;
        R_P_max = P_max;
        R_Q_min = Q_max*(V_k - (1-del_perm))/(V_min - (1-del_perm));
        R_Q_max = Q_max;     
    elseif V_k >= 1- del_perm && V_k <= 1+ del_perm
        R_P_min = P_min;
        R_P_max = P_max;
        R_Q_min = Q_min;
        R_Q_max = Q_max;
    elseif V_k > 1+ del_perm && V_k <= V_max
        R_P_min = P_min;
        R_P_max = P_max;
        R_Q_min = Q_min;
        R_Q_max = Q_min*(V_k - (1+del_perm))/(V_max - (1+del_perm));
    elseif V_k > V_max
        R_P_min = P_min;
        R_P_max = P_max;
        R_Q_min = Q_min;
        R_Q_max = Q_min;
    end
elseif mode == 5        %Q(V) by ANRC
    if V_k < V_min
        R_P_min = P_min;
        R_P_max = P_max;
        R_Q_min = 0;
        R_Q_max = Q_max;
    elseif V_k >= V_min && V_k < 1- del_perm
        R_P_min = P_min;
        R_P_max = P_max;
        R_Q_min = Q_min*(V_min - V_k)/(V_min - (1-del_perm));
        R_Q_max = Q_max;     
    elseif V_k >= 1- del_perm && V_k <= 1+ del_perm
        R_P_min = P_min;
        R_P_max = P_max;
        R_Q_min = Q_min;
        R_Q_max = Q_max;
    elseif V_k > 1+ del_perm && V_k <= V_max
        R_P_min = P_min;
        R_P_max = P_max;
        R_Q_min = Q_min;
        R_Q_max = Q_max*(V_max - V_k)/(V_max - (1+del_perm));
    elseif V_k > V_max
        R_P_min = P_min;
        R_P_max = P_max;
        R_Q_min = Q_min;
        R_Q_max = 0;
    end
    
else
    R_P_min = P_min;
    R_P_max = P_max;
    R_Q_min = Q_min;
    R_Q_max = Q_max;
end

% Prior to using this function fix P_min, P_max, Q_min, Q_max